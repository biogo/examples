// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package igor

import (
	"log"

	"github.com/biogo/biogo/align/pals"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

// store is a pals.Contig internment implementation.
type store map[pals.Contig]pals.Contig

// intern returns an interned copy of the parameter.
func (is store) intern(c pals.Contig) pals.Contig {
	if c == "" {
		return ""
	}
	t, ok := is[c]
	if ok {
		return t
	}
	is[c] = c
	return c
}

// Piles reads the features in the input gff.Reader and applies pals.Piler analysis
// using the specified overlap and pair filter function. The features in the input
// must satisfy pals.ExpandFeature restrictions.
func Piles(in *gff.Reader, overlap int, pf pals.PairFilter, l *log.Logger, freq int) ([]*pals.Pile, error) {
	piler := pals.NewPiler(overlap)
	contigs := make(store)

	var n int

	sc := featio.NewScanner(in)
	for sc.Next() {
		rep := sc.Feat().(*gff.Feature)

		p, err := pals.ExpandFeature(rep)
		if err != nil {
			return nil, err
		}
		p.A.Loc = contigs.intern(p.A.Loc.(pals.Contig))
		p.B.Loc = contigs.intern(p.B.Loc.(pals.Contig))

		piler.Add(p)
		n++
		if l != nil && freq != 0 && n%freq == 0 {
			l.Printf("added %d images", n)
		}
	}
	if err := sc.Error(); err != nil {
		return nil, err
	}

	piler.Logger = l
	piler.LogFreq = freq
	return piler.Piles(pf), nil
}
