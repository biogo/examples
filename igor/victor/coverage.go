// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"os"
	"sort"
	"strings"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/store/interval"
	"github.com/biogo/store/step"
)

type feature struct {
	Chr    string     `json:"C"`
	Start  int        `json:"S"`
	End    int        `json:"E"`
	Orient seq.Strand `json:"O"`

	id uintptr
}

func (f feature) Overlap(b interval.IntRange) bool {
	return f.End >= b.Start && f.Start <= b.End
}
func (f feature) Range() interval.IntRange {
	return interval.IntRange{Start: f.Start, End: f.End}
}
func (f feature) ID() uintptr { return f.id }

// family is a collection of features that share sequence identity.
type family struct {
	id      int64     // node ID for graph analysis
	members []feature // features within the family
	length  int       // total genome coverage by members of the family
}

// byMembers is a sort helper that sorts a []family descending by the number of
// members in the families.
type byMembers []family

func (f byMembers) Len() int           { return len(f) }
func (f byMembers) Less(i, j int) bool { return len(f[i].members) > len(f[j].members) }
func (f byMembers) Swap(i, j int)      { f[i], f[j] = f[j], f[i] }

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

// pair is a [2]bool type satisfying the step.Equaler interface.
type pair [2]bool

// Equal returns whether p equals e. Equal assumes the underlying type of e is pair.
func (p pair) Equal(e step.Equaler) bool {
	return p == e.(pair)
}

// intersection returns the fractions of intersection of a and b, where
// upper is |a∩b|/min(|a|,|b|) and lower is |a∩b|/max(|a|,|b|).
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

// flattenFamily returns an interval forest representing the location
// of features in the family such that intervals are disjoint.
func flattenFamily(f family) (map[string]*interval.IntTree, error) {
	if len(f.members) == 0 {
		return nil, nil
	}

	trees := make(map[string]*interval.IntTree)
	var id uintptr
	for _, m := range f.members {
		t, ok := trees[m.Chr]
		if !ok {
			t = &interval.IntTree{}
			trees[m.Chr] = t
		}
		m.id = id
		add := m
		var del []interval.IntInterface
		hits := t.Get(m)
		for _, h := range hits {
			hr := h.Range()
			if m.Start < hr.Start || hr.End < m.End {
				del = append(del, h)
			}
			if hr.Start < add.Start {
				add.Start = hr.Start
			}
			if add.End < hr.End {
				add.End = hr.End
			}
		}
		for _, d := range del {
			t.Delete(d, true)
		}
		if len(hits) == 0 || len(del) != 0 {
			err := t.Insert(add, true)
			if err != nil {
				return nil, err
			}
		}
		t.AdjustRanges()
		id++
	}

	return trees, nil
}

type nameSupport struct {
	name     string
	coverage float64
}

type nameSupports []nameSupport

func (n nameSupports) String() string {
	var buf strings.Builder
	for i, v := range n {
		if i != 0 {
			buf.WriteByte(',')
		}
		fmt.Fprintf(&buf, "%q=%f", v.name, v.coverage)
	}
	return buf.String()
}

// nameCoverage returns a slice of name support for for the family coverage
// obtained from flattenFamily in famcov using annotation in annots.
func nameCoverage(famcov, annots map[string]*interval.IntTree, normalise bool) []nameSupport {
	names := make(map[string]float64)
	var size float64
	if normalise {
		for _, intervals := range famcov {
			intervals.Do(func(iv interval.IntInterface) (done bool) {
				r := iv.Range()
				size += float64(r.End - r.Start)
				return
			})
		}
	} else {
		size = 1
	}
	for chr, intervals := range famcov {
		ann, ok := annots[chr]
		if !ok {
			continue
		}

		intervals.Do(func(q interval.IntInterface) (done bool) {
			qr := q.Range()
			ann.DoMatching(func(a interval.IntInterface) (done bool) {
				ar := a.Range()

				n := a.(annotation).FeatAttributes.Get("Repeat")
				first := false
				for i, r := range n {
					if r == ' ' {
						if first {
							n = n[:i]
							break
						}
						first = true
					}
				}

				names[n] += float64(min(qr.End, ar.End) - max(qr.Start, ar.Start))
				return
			}, q)
			return
		})
	}

	ns := make([]nameSupport, 0, len(names))
	for n, c := range names {
		ns = append(ns, nameSupport{name: n, coverage: c / size})
	}
	sort.Sort(bySupport(ns))

	return ns
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

// bySupport is a sort helper that sorts a []nameSupport descending by the base coverage
// of the name in the family.
type bySupport []nameSupport

func (n bySupport) Len() int           { return len(n) }
func (n bySupport) Less(i, j int) bool { return n[i].coverage > n[j].coverage }
func (n bySupport) Swap(i, j int)      { n[i], n[j] = n[j], n[i] }

// annotation is an gff.Feature that satisfies interval.IntInterface.
type annotation struct {
	*gff.Feature
	id uintptr
}

func (a annotation) ID() uintptr { return a.id }

func (a annotation) Overlap(b interval.IntRange) bool {
	return a.FeatEnd >= b.Start && a.FeatStart <= b.End
}

func (a annotation) Range() interval.IntRange {
	return interval.IntRange{Start: a.FeatStart, End: a.FeatEnd}
}

// readAnnotation reads in a BE RepeatMasker GFF and
// returns an interval forest of annotations.
func readAnnotations(path string) (map[string]*interval.IntTree, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	annot := make(map[string]*interval.IntTree)
	sc := featio.NewScanner(gff.NewReader(f))
	var id uintptr
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		t, ok := annot[f.SeqName]
		if !ok {
			t = &interval.IntTree{}
			annot[f.SeqName] = t
		}
		err := t.Insert(annotation{Feature: f, id: id}, true)
		if err != nil {
			return nil, err
		}
		id++
	}
	err = sc.Error()
	if err != nil {
		return nil, err
	}
	for _, t := range annot {
		t.AdjustRanges()
	}
	return annot, nil
}
