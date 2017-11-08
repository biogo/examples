// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// victor is a post processor for grouping families defined by igor.
package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime"
	"sort"
	"sync"

	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/store/step"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/community"
	"gonum.org/v1/gonum/graph/encoding"
	"gonum.org/v1/gonum/graph/encoding/dot"
	"gonum.org/v1/gonum/graph/network"
	"gonum.org/v1/gonum/graph/simple"
	"gonum.org/v1/gonum/graph/topo"
)

var (
	in         = flag.String("in", "", "Specifies the input json file name.")
	dotOut     = flag.String("dot", "", "Specifies the output DOT file name.")
	thresh     = flag.Float64("thresh", 0.05, "Specifies minimum family intersection to report.")
	resolution = flag.Float64("resolution", 1, "Specifies the resolution for cluster modularisation.")
	minFam     = flag.Int("min", 0, "Specify the minimum number of members in a family to include (if 0 no limit).")
	cliques    = flag.Bool("cliques", false, "Find cliques in non-clique clusters.")
	threads    = flag.Int("threads", 0, "Specify the number of parallel connection threads (if 0 use GOMAXPROCS).")
)

func main() {
	flag.Parse()
	if *in == "" {
		flag.Usage()
		os.Exit(0)
	}
	if *threads == 0 {
		*threads = runtime.GOMAXPROCS(0)
	}

	f, err := os.Open(*in)
	if err != nil {
		log.Fatalf("failed reading %q: %v", *in, err)
	}
	defer f.Close()
	r := bufio.NewReader(f)

	var families []family
	for i := 0; ; i++ {
		l, err := r.ReadBytes('\n')
		if err != nil {
			break
		}
		var v []feature
		err = json.Unmarshal(l, &v)
		if err != nil {
			log.Fatalf("failed unmarshaling json for family %d: %v", i, err)
		}
		if *minFam != 0 && len(v) < *minFam {
			continue
		}
		fam := family{id: int64(i), members: v, length: length(v)}

		families = append(families, fam)
	}
	sort.Sort(byMembers(families))

	c := connector{limit: make(chan struct{}, *threads)}
	edges := c.edgesFor(families, *thresh)

	const minSubClique = 3
	grps := groups(families, edges, *resolution, minSubClique, *cliques)

	clusterIdentity := make(map[int64]int64)
	cliqueIdentity := make(map[int64][]int64)
	cliqueMemberships := make(map[int64]int64)

	for _, g := range grps {
		// Collate counts for clique memberships. We cannot do
		// this one group at a time; a member of a group can be
		// a clique member of another group since they are in
		// potentially in connection with other groups.
		for _, clique := range g.cliques {
			for _, m := range clique {
				cliqueMemberships[m]++
			}
		}
	}
	for _, g := range grps {
		fmt.Fprintf(os.Stderr, "clique=%t", g.isClique)
		for _, m := range g.members {
			fmt.Fprintf(os.Stderr, " %d", m.id)
			clusterIdentity[m.id] = g.pageRank[0].id
			if g.isClique {
				cliqueMemberships[m.id]++
				cliqueIdentity[m.id] = []int64{g.pageRank[0].id}
			}
		}
		if len(g.cliques) != 0 {
			fmt.Fprintf(os.Stderr, " (%d+)-cliquesIn=%v", minSubClique, g.cliques)
		}
		for _, clique := range g.cliques {
			// Make PageRanked version of clique.
			cliqueHas := make(map[int64]bool)
			for _, m := range clique {
				cliqueHas[m] = true
			}
			clique = make([]int64, 0, len(clique))
			for _, m := range g.pageRank {
				if cliqueHas[m.id] {
					clique = append(clique, m.id)
				}
			}

			// Annotate families as meaningfully but concisely as possible.
			unique := cliqueMemberships[clique[0]] == 1
			for i, m := range clique {
				if cliqueMemberships[m] == 1 {
					if unique {
						cliqueIdentity[m] = clique[:1]
					} else {
						cliqueIdentity[m] = clique
					}
				} else {
					cliqueIdentity[m] = clique[i : i+1]
				}
			}
		}
		fmt.Fprintf(os.Stderr, " PageRank=%+v\n", g.pageRank)
	}
	for i, e := range edges {
		if clustID, isClustered := clusterIdentity[e.from.id]; isClustered {
			edges[i].from.cluster = clustID
		}
		if clustID, isClustered := clusterIdentity[e.to.id]; isClustered {
			edges[i].to.cluster = clustID
		}
	}
	if *dotOut != "" {
		writeDOT(*dotOut, edges)
	}

	b := bufio.NewWriter(os.Stdout)
	defer b.Flush()
	w := gff.NewWriter(b, 60, false)
	ft := &gff.Feature{
		Source:  "igor/victor",
		Feature: "repeat",
		FeatAttributes: gff.Attributes{
			{Tag: "Family"},
			{Tag: "Cluster"},
			{Tag: "Clique"},
		},
	}
	for _, fam := range families {
		clustID, isClustered := clusterIdentity[fam.id]
		for _, m := range fam.members {
			ft.SeqName = m.Chr
			ft.FeatStart = m.Start
			ft.FeatEnd = m.End
			ft.FeatStrand = m.Orient
			ft.FeatFrame = gff.NoFrame
			ft.FeatAttributes[0].Value = fmt.Sprint(fam.id)
			if isClustered {
				ft.FeatAttributes = ft.FeatAttributes[:3]
				ft.FeatAttributes[1].Value = fmt.Sprint(clustID)
				switch id := cliqueIdentity[fam.id]; {
				case id == nil:
					ft.FeatAttributes = ft.FeatAttributes[:2]
				case cliqueMemberships[fam.id] == 1:
					ft.FeatAttributes[2].Value = dotted(id)
				default:
					ft.FeatAttributes[2].Value = fmt.Sprintf("%d*", id[0])
				}
			} else {
				ft.FeatAttributes = ft.FeatAttributes[:1]
			}
			_, err := w.Write(ft)
			if err != nil {
				log.Fatalf("error: %v", err)
			}
		}
	}
}

type feature struct {
	Chr    string     `json:"C"`
	Start  int        `json:"S"`
	End    int        `json:"E"`
	Orient seq.Strand `json:"O"`
}

type family struct {
	id      int64
	members []feature
	length  int
}

type byMembers []family

func (f byMembers) Len() int           { return len(f) }
func (f byMembers) Less(i, j int) bool { return len(f[i].members) > len(f[j].members) }
func (f byMembers) Swap(i, j int)      { f[i], f[j] = f[j], f[i] }

type node struct {
	id      int64
	cluster int64
	members int
}

var _ encoding.Attributer = node{}

func (n node) ID() int64 { return n.id }
func (n node) Attributes() []encoding.Attribute {
	if n.cluster == -1 {
		return []encoding.Attribute{{"members", fmt.Sprint(n.members)}}
	}
	return []encoding.Attribute{
		{"cluster", fmt.Sprint(n.cluster)},
		{"members", fmt.Sprint(n.members)},
	}
}

type edge struct {
	from, to node
	weight   float64
}

var _ encoding.Attributer = edge{}

func (e edge) From() graph.Node { return e.from }
func (e edge) To() graph.Node   { return e.to }
func (e edge) Weight() float64  { return e.weight }
func (e edge) Attributes() []encoding.Attribute {
	return []encoding.Attribute{{"weight", fmt.Sprint(e.weight)}}
}

// stepBool is a bool type satisfying the step.Equaler interface.
type stepBool bool

// Equal returns whether b equals e. Equal assumes the underlying type of e is a stepBool.
func (b stepBool) Equal(e step.Equaler) bool {
	return b == e.(stepBool)
}

// length returns the number of covered bases in v.
func length(v []feature) int {
	vecs := make(map[string]*step.Vector)
	for _, f := range v {
		vec, ok := vecs[f.Chr]
		if !ok {
			var err error
			vec, err = step.New(f.Start, f.End, stepBool(false))
			if err != nil {
				panic(err)
			}
			vec.Relaxed = true
			vecs[f.Chr] = vec
		}
		vec.SetRange(f.Start, f.End, stepBool(true))
	}
	var len int
	for _, vec := range vecs {
		vec.Do(func(start, end int, e step.Equaler) {
			if e.(stepBool) {
				len += end - start
			}
		})
	}
	return len
}

// connector handles parallel analysis of family intersections.
type connector struct {
	wg sync.WaitGroup

	mu    sync.Mutex
	edges []edge

	// limit specifies the maximum number
	// of concurrent intersection calls.
	limit chan struct{}
}

// acquire gets an available worker thread.
func (c *connector) acquire() {
	c.wg.Add(1)
	c.limit <- struct{}{}
}

// release puts pack a worker thread.
func (c *connector) release() {
	<-c.limit
	c.wg.Done()
}

// connect adds e to the store of edges.
func (c *connector) connect(e edge) {
	c.mu.Lock()
	fmt.Fprintln(os.Stderr, e.from.id, e.to.id, e.weight)
	c.edges = append(c.edges, e)
	c.mu.Unlock()
}

// edgesFor returns the edges that exist between families in f where
// the intersection is greater than or equal to thresh.
func (c *connector) edgesFor(f []family, thresh float64) []edge {
	for i, a := range f[:len(f)-1] {
		for _, b := range f[i+1:] {
			a := a
			b := b
			c.acquire()
			go func() {
				defer c.release()
				upper, lower := intersection(a, b)
				if upper < thresh {
					return
				}

				// Edges indicate connection from the shorter
				// family to the longer family, so ensure this
				// is the state now.
				if a.length > b.length {
					a, b = b, a
				}

				c.connect(edge{
					from:   node{id: a.id, cluster: -1, members: len(a.members)},
					to:     node{id: b.id, cluster: -1, members: len(b.members)},
					weight: upper,
				})

				if lower < thresh {
					return
				}

				c.connect(edge{
					from:   node{id: b.id, cluster: -1, members: len(b.members)},
					to:     node{id: a.id, cluster: -1, members: len(a.members)},
					weight: lower,
				})
			}()
		}
	}
	c.wg.Wait()

	c.mu.Lock()
	defer c.mu.Unlock()

	return c.edges
}

// pair is a [2]bool type satisfying the step.Equaler interface.
type pair [2]bool

// Equal returns whether p equals e. Equal assumes the underlying type of e is pair.
func (p pair) Equal(e step.Equaler) bool {
	return p == e.(pair)
}

func intersection(a, b family) (upper, lower float64) {
	// TODO(kortschak): Consider orientation agreement.
	vecs := make(map[string]*step.Vector)
	for i, v := range []family{a, b} {
		for _, f := range v.members {
			vec, ok := vecs[f.Chr]
			if !ok {
				var err error
				vec, err = step.New(f.Start, f.End, pair{})
				if err != nil {
					panic(err)
				}
				vec.Relaxed = true
				vecs[f.Chr] = vec
			}
			err := vec.ApplyRange(f.Start, f.End, func(e step.Equaler) step.Equaler {
				p := e.(pair)
				p[i] = true
				return p
			})
			if err != nil {
				panic(err)
			}
		}
	}
	var (
		aLen, bLen int
		intersect  int
	)
	for _, vec := range vecs {
		vec.Do(func(start, end int, e step.Equaler) {
			p := e.(pair)
			if p[0] {
				aLen += end - start
			}
			if p[1] {
				bLen += end - start
			}
			if p[0] && p[1] {
				intersect += end - start
			}
		})
	}
	if aLen != a.length || bLen != b.length {
		panic("length mismatch")
	}

	upper = float64(intersect) / math.Min(float64(a.length), float64(b.length))
	lower = float64(intersect) / math.Max(float64(a.length), float64(b.length))
	return upper, lower
}

func dotted(id []int64) string {
	var buf bytes.Buffer
	for i, e := range id {
		if i != 0 {
			fmt.Fprint(&buf, ".")
		}
		fmt.Fprint(&buf, e)
	}
	return buf.String()
}

func writeDOT(file string, edges []edge) {
	g := simple.NewWeightedDirectedGraph(0, math.Inf(1))
	for _, e := range edges {
		for _, n := range []graph.Node{e.From(), e.To()} {
			if !g.Has(n) {
				g.AddNode(n)
			}
		}
		g.SetWeightedEdge(e)
	}

	f, err := os.Create(*dotOut)
	if err != nil {
		log.Printf("failed to create %q DOT output file: %v", *dotOut, err)
		return
	}
	defer f.Close()
	b, err := dot.Marshal(g, "", "", "  ", false)
	if err != nil {
		log.Printf("failed to create DOT bytes: %v", err)
		return
	}
	_, err = f.Write(b)
	if err != nil {
		log.Printf("failed to write DOT: %v", err)
	}
}

type group struct {
	members  []family
	isClique bool
	cliques  [][]int64
	pageRank ranks
}

func groups(fams []family, edges []edge, resolution float64, minSubClique int, cliques bool) []group {
	g := simple.NewWeightedDirectedGraph(0, 0)
	for _, e := range edges {
		for _, n := range []graph.Node{e.From(), e.To()} {
			if !g.Has(n) {
				g.AddNode(n)
			}
		}
		g.SetWeightedEdge(e)
	}

	familyIndexOf := make(map[int64]int, len(fams))
	for i, f := range fams {
		familyIndexOf[f.id] = i
	}
	var grps []group
	r := community.Modularize(graph.Undirect{G: g}, resolution, rand.New(rand.NewSource(1)))
	for _, c := range r.Communities() {
		var grp group
		for _, n := range c {
			grp.members = append(grp.members, fams[familyIndexOf[n.ID()]])
		}
		if len(grp.members) == 2 || edgesIn(g, c)*2 == len(c)*(len(c)-1) {
			grp.isClique = true
		} else if cliques {
			grp.cliques = cliquesIn(grp, edges, minSubClique)
		}
		if len(grp.members) > 1 {
			grp.pageRank = ranksOf(grp, edges)
		}

		grps = append(grps, grp)
	}

	return grps
}

type intset map[int64]struct{}

func (s intset) add(i int64) {
	s[i] = struct{}{}
}

func (s intset) has(i int64) bool {
	_, ok := s[i]
	return ok
}

type twoset map[[2]int64]struct{}

func (s twoset) add(i, j int64) {
	if i > j {
		i, j = j, i
	}
	s[[2]int64{i, j}] = struct{}{}
}

func edgesIn(g graph.Directed, n []graph.Node) int {
	in := make(intset)
	for _, u := range n {
		in.add(u.ID())
	}
	seen := make(twoset)
	// We could use graph.Undirect here, but the
	// overhead increases and we don't actually
	// need all the nodes, just the edges.
	for _, u := range n {
		uid := u.ID()
		for _, v := range g.From(u) {
			vid := v.ID()
			if !in.has(vid) {
				continue
			}
			seen.add(uid, vid)
		}
		for _, v := range g.To(u) {
			vid := v.ID()
			if !in.has(vid) {
				continue
			}
			seen.add(uid, vid)
		}
	}
	return len(seen)
}

func cliquesIn(grp group, edges []edge, min int) [][]int64 {
	members := make(intset)
	for _, fam := range grp.members {
		members.add(fam.id)
	}

	g := simple.NewUndirectedGraph()
outer:
	for _, e := range edges {
		for _, n := range []graph.Node{e.From(), e.To()} {
			if !members.has(n.ID()) {
				continue outer
			}
		}
		for _, n := range []graph.Node{e.From(), e.To()} {
			if !g.Has(n) {
				g.AddNode(n)
			}
		}
		g.SetEdge(e)
	}

	clqs := topo.BronKerbosch(g)
	var cliqueIDs [][]int64
	for _, clq := range clqs {
		if len(clq) < min {
			continue
		}
		ids := make([]int64, 0, len(clq))
		for _, n := range clq {
			ids = append(ids, n.ID())
		}
		cliqueIDs = append(cliqueIDs, ids)
	}

	return cliqueIDs
}

func ranksOf(grp group, edges []edge) ranks {
	members := make(intset)
	for _, fam := range grp.members {
		members.add(fam.id)
	}

	g := simple.NewWeightedDirectedGraph(0, math.Inf(1))
outer:
	for _, e := range edges {
		for _, n := range []graph.Node{e.From(), e.To()} {
			if !members.has(n.ID()) {
				continue outer
			}
		}
		for _, n := range []graph.Node{e.From(), e.To()} {
			if !g.Has(n) {
				g.AddNode(n)
			}
		}
		g.SetWeightedEdge(e)
	}

	r := network.PageRank(g, 0.85, 1e-6)
	o := make(ranks, 0, len(r))
	for id, rnk := range r {
		o = append(o, rank{id: id, rank: rnk})
	}
	sort.Sort(o)
	return o
}

type rank struct {
	id   int64
	rank float64
}

type ranks []rank

func (o ranks) Len() int { return len(o) }
func (o ranks) Less(i, j int) bool {
	return o[i].rank > o[j].rank || (o[i].rank == o[j].rank && o[i].id < o[j].id)
}
func (o ranks) Swap(i, j int) { o[i], o[j] = o[j], o[i] }
