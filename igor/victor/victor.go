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
	"os"
	"runtime"
	"sort"

	"github.com/biogo/biogo/io/featio/gff"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/encoding/dot"
	"gonum.org/v1/gonum/graph/simple"
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

// writeDot writes the graph in edges to the named file in DOT.
func writeDOT(file string, edges []edge) {
	g := simple.NewWeightedDirectedGraph(0, math.Inf(1))
	for _, e := range edges {
		for _, n := range []graph.Node{e.From(), e.To()} {
			if !g.Has(n.ID()) {
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

// dotted formats the IDs in id separated by '.'.
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
