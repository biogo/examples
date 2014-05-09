// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bufio"
	"code.google.com/p/biogo.graph"
	"code.google.com/p/biogo.store/interval"
	"code.google.com/p/biogo/seq"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"unsafe"
)

type Trees struct {
	graph.Node
	Intervals map[string]*interval.IntTree
}

func NewTrees(v []*feat) *Trees {
	ts := &Trees{Intervals: make(map[string]*interval.IntTree)}
	for _, i := range v {
		ts.Insert(i, i.Chr, true)
	}
	ts.AdjustRanges()
	return ts
}

func (ts *Trees) Len() int {
	var l int
	for _, t := range ts.Intervals {
		l += t.Len()
	}
	return l
}

func (ts *Trees) Segments() []string {
	s := make([]string, 0, len(ts.Intervals))
	for k := range ts.Intervals {
		s = append(s, k)
	}
	return s
}

func (ts *Trees) Insert(i interval.IntInterface, seg string, fast bool) error {
	t, ok := ts.Intervals[seg]
	if !ok {
		t = &interval.IntTree{}
		ts.Intervals[seg] = t
	}
	return t.Insert(i, fast)
}

func (ts *Trees) Delete(e interval.IntInterface, seg string, fast bool) error {
	t, ok := ts.Intervals[seg]
	if !ok {
		return nil
	}
	return t.Delete(e, fast)
}

func (ts *Trees) Get(q interval.IntInterface, seg string) []interval.IntInterface {
	t, ok := ts.Intervals[seg]
	if !ok {
		return nil
	}
	return t.Get(q)
}

func (ts *Trees) DoMatching(fn interval.IntOperation, q interval.IntInterface, seg string) bool {
	t, ok := ts.Intervals[seg]
	if !ok {
		return false
	}
	return t.DoMatching(fn, q)
}

func (ts *Trees) Do(fn interval.IntOperation, seg string) bool {
	t, ok := ts.Intervals[seg]
	if !ok {
		return false
	}
	return t.Do(fn)
}

func (ts *Trees) AdjustRanges() {
	for _, t := range ts.Intervals {
		t.AdjustRanges()
	}
}

type feat struct {
	Chr    string     `json:"C"`
	Start  int        `json:"S"`
	End    int        `json:"E"`
	Orient seq.Strand `json:"O"`
}

func (f *feat) Overlap(b interval.IntRange) bool { return f.End > b.Start && f.Start < b.End }
func (f *feat) ID() uintptr                      { return uintptr(unsafe.Pointer(f)) }
func (f *feat) Range() interval.IntRange         { return interval.IntRange{f.Start, f.End} }

type strandEdge struct {
	graph.Edge
	Strand seq.Strand
}

var (
	maxFam  int
	epsilon float64
)

func main() {
	flag.IntVar(&maxFam, "maxFam", 0, "maxFam indicates maximum family size considered (0 == no limit).")
	flag.Float64Var(&epsilon, "epsilon", 0.0225, "Tolerance for clustering.")
	flag.Parse()

	if len(flag.Args()) < 1 {
		fmt.Fprintln(os.Stderr, "Need input file.")
		os.Exit(1)
	}
	if len(flag.Args()) < 2 {
		fmt.Fprintln(os.Stderr, "Nothing to do.")
		os.Exit(0)
	}

	var (
		g   = graph.NewUndirected()
		bad int
	)
	{
		var tss []*Trees
		for _, n := range flag.Args() {
			f, err := os.Open(n)
			if err != nil {
				fmt.Fprintf(os.Stdout, "Error: %v\n", err)
				os.Exit(1)
			}
			b := bufio.NewReader(f)

			var (
				tas []*Trees
				ori = make(map[struct{ k, j int }]seq.Strand)
			)
			for j := 0; ; j++ {
				l, err := b.ReadBytes('\n')
				if err != nil {
					break
				}
				v := []*feat{}
				err = json.Unmarshal(l, &v)
				if err != nil {
					fmt.Fprintf(os.Stdout, "Error: %v\n", err)
					os.Exit(1)
				}
				if maxFam != 0 && len(v) > maxFam {
					continue
				}

				jn := NewTrees(v)
				jn.Node = g.NewNode()
				g.Add(jn)
				tas = append(tas, jn)

				if tss != nil {
					// Search tss for good matches with the current family...
					for _, i := range v {
						for k, ts := range tss {
							ts.DoMatching(func(iv interval.IntInterface) (done bool) {
								p := iv.(*feat)
								if isClose(p, i, epsilon) {
									o, ok := ori[struct{ k, j int }{k, j}]
									if !ok {
										ori[struct{ k, j int }{k, j}] = p.Orient * i.Orient
									} else if o != p.Orient*i.Orient {
										bad++
										fmt.Fprintln(os.Stderr, "#### BAD ORIENTATION ####")
									}

									con, err := g.Connected(ts, jn)
									if err != nil {
										panic(err)
									}
									if !con {
										g.ConnectWith(ts, jn, strandEdge{Edge: graph.NewEdge(), Strand: p.Orient * i.Orient})
									}
								}
								return
							}, i, i.Chr)
						}
					}
				}
			}
			f.Close()
			if tss == nil {
				tss = tas
			} else {
				tss = append(tss, tas...)
			}
		}
	}

	cc := g.ConnectedComponents(graph.EdgeFilter(func(e graph.Edge) bool {
		// We need to correct orientation here.
		return true
	}))
	fmt.Printf("Bad orientation connections: %d G=%v Connected components: %d\n", bad, g, len(cc))

	for i, c := range cc {
		ts := c[0].(*Trees)
		for j, fi := range c[1:] {
			cfi := fi.(*Trees)
			for _, s := range cfi.Segments() {
				cfi.Do(func(e interval.IntInterface) (done bool) {
					var exact bool
					for _, m := range ts.Get(e, s) {
						if m.Range() == e.Range() {
							exact = true
						}
					}
					if !exact {
						ts.Insert(e, s, true)
					}
					return
				}, s)
			}
			c[j+1] = nil
		}
		cc[i] = cc[i][:1]
		fmt.Printf("Component %d: %d\n", i, ts.Len())
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func isClose(a, b *feat, thresh float64) bool {
	s, e := min(a.Start, b.Start), max(a.End, b.End)
	l := float64(e-s) / 4
	cut := thresh * l * l
	ds, de := float64(a.Start-b.Start), float64(a.End-b.End)
	return ds*ds+de*de < cut
}
