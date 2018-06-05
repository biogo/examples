// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math"

	"github.com/biogo/biogo/seq"
	"github.com/biogo/store/step"
)

type feature struct {
	Chr    string     `json:"C"`
	Start  int        `json:"S"`
	End    int        `json:"E"`
	Orient seq.Strand `json:"O"`
}

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

// interstion returns the fractions of intersection of a and b, where
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
