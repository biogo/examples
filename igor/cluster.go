package main

import (
	"code.google.com/p/biogo.cluster/kmeans"
	"code.google.com/p/biogo/align/pals"
	"code.google.com/p/biogo/feat"
	"fmt"
)

type Pairs []*pals.Pair

func (p Pairs) Len() int                    { return len(p) }
func (p Pairs) Values(i int) (x, y float64) { return float64(p[i].A.Start()), float64(p[i].A.End()) }

type Cluster struct {
	From  int
	To    int
	Loc   feat.Feature
	Piles []*pals.Pile
}

func (cl *Cluster) Name() string {
	return fmt.Sprintf("%s[%d,%d)", cl.Loc.Name(), cl.From, cl.To)
}

// Description returns the string "pals feature".
func (cl *Cluster) Description() string    { return "cluster" }
func (cl *Cluster) Start() int             { return cl.From }
func (cl *Cluster) End() int               { return cl.To }
func (cl *Cluster) Len() int               { return cl.To - cl.From }
func (cl *Cluster) Location() feat.Feature { return cl.Loc }

func (cl *Cluster) String() string {
	return fmt.Sprintf("{%s[%d,%d): %v}", cl.Loc.Name(), cl.From, cl.To, cl.Piles)
}

func BuildCluster(p *pals.Pile, epsilon float64, effort int) *Cluster {
	km := kmeans.NewKmeans(Pairs(p.Images))

	values := km.Values()
	cut := make([]float64, len(values))
	for i, v := range values {
		l := epsilon * (v.Y() - v.X())
		cut[i] = l * l
	}

	for k := 1; k <= len(p.Images); k++ {
		var e int
		if k == 1 {
			e = 1
		} else {
			e = effort
		}
	L:
		for a := 0; a < e; a++ {
			km.Seed(k)
			km.Cluster()
			centers := km.Means()
			for i, v := range values {
				dx, dy := centers[v.Cluster()].X()-v.X(), centers[v.Cluster()].Y()-v.Y()
				ok := dx*dx+dy*dy < cut[i]
				if !ok {
					continue L
				}
			}

			var (
				pc = make([]*pals.Pile, k)
				cl = &Cluster{
					From:  p.Start(),
					To:    p.End(),
					Loc:   p.Location(),
					Piles: pc,
				}
			)
			for ci, c := range km.Clusters() {
				pc[ci] = &pals.Pile{
					From:   int(centers[ci].X()),
					To:     int(centers[ci].Y()),
					Loc:    cl,
					Images: make([]*pals.Pair, centers[ci].Count()),
				}
				for pi, i := range c {
					im := p.Images[i]
					im.A.(*pals.Feature).Loc = pc[ci]
					pc[ci].Images[pi] = im
				}
			}

			return cl
		}
	}

	panic("cannot reach")
}
