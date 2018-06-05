// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"os"
	"sort"
	"sync"

	"golang.org/x/exp/rand"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/community"
	"gonum.org/v1/gonum/graph/encoding"
	"gonum.org/v1/gonum/graph/network"
	"gonum.org/v1/gonum/graph/simple"
	"gonum.org/v1/gonum/graph/topo"
)

// node is a graph node with DOT support.
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

// edge is a graph edge.
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

// group is a collection of families related by shared genomic
// locations
type group struct {
	// members is the collection of families consituting the group.
	members []family
	// isClique indicates that the group is fully connected.
	isClique bool
	// cliques is the node IDs of subcliques within the group.
	cliques [][]int64
	// page rank is the rank of importance of each family.
	// Importance here is a measure of how much each family
	// covers relative to its overlapping partners.
	pageRank ranks
}

// groups returns a collection of grouped families.
//
// The families are grouped by constructing a weighted graph of family
// relationships defined by total family interval overlap in edges.
// Clustering is performed using Louvain community detection at the
// specified resolution.
//
// If cliques is true, cliques at least minSubClique within each group
// will be identified.
func groups(fams []family, edges []edge, resolution float64, minSubClique int, cliques bool) []group {
	g := simple.NewWeightedDirectedGraph(0, 0)
	for _, e := range edges {
		for _, n := range []graph.Node{e.From(), e.To()} {
			if !g.Has(n.ID()) {
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

// inset is a simple integer set.
type intset map[int64]struct{}

func (s intset) add(i int64) {
	s[i] = struct{}{}
}

func (s intset) has(i int64) bool {
	_, ok := s[i]
	return ok
}

// twoset is a simples tuple set.
type twoset map[[2]int64]struct{}

func (s twoset) add(i, j int64) {
	if i > j {
		i, j = j, i
	}
	s[[2]int64{i, j}] = struct{}{}
}

// edgesIn returns the number of edges in the graph g vertex-induced with n.
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
		for _, v := range g.From(uid) {
			vid := v.ID()
			if !in.has(vid) {
				continue
			}
			seen.add(uid, vid)
		}
		for _, v := range g.To(uid) {
			vid := v.ID()
			if !in.has(vid) {
				continue
			}
			seen.add(uid, vid)
		}
	}
	return len(seen)
}

// cliquesIn returns the sub-cliques with at least min nodes
// present in the graph defined by edges vertex induced by
// families within grp.
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
			if !g.Has(n.ID()) {
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

// ranksOf returns the Page Rank ranking of families in group
// where
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
			if !g.Has(n.ID()) {
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

// rank is an element in a ordered key value store.
type rank struct {
	id   int64
	rank float64
}

// ranks is a collection of ranks that sorts descending by rank, or
// ascending by id when ranks are equal.
type ranks []rank

func (o ranks) Len() int { return len(o) }
func (o ranks) Less(i, j int) bool {
	return o[i].rank > o[j].rank || (o[i].rank == o[j].rank && o[i].id < o[j].id)
}
func (o ranks) Swap(i, j int) { o[i], o[j] = o[j], o[i] }
