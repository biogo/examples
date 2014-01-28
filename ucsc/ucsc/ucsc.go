// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package ucsc provides a linear sequence type that parses UCSC header data into
// the sequence annotation data.
package ucsc

import (
	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/feat"
	"code.google.com/p/biogo/seq"
	"code.google.com/p/biogo/seq/linear"

	"strconv"
	"strings"
)

// A Chr is a feat.Feature marking a chromosome.
type Chr string

func (c Chr) Start() int             { return 0 }
func (c Chr) End() int               { return 0 }
func (c Chr) Len() int               { return 0 }
func (c Chr) Name() string           { return string(c) }
func (c Chr) Description() string    { return "chromosome" }
func (c Chr) Location() feat.Feature { return nil }

// Seq modifies the behaviour of linear.Seq so that the description is parsed
// according to the UCSC format.
type Seq struct {
	*linear.Seq
}

// NewSeq returns a new Seq.
func NewSeq(id string, b []alphabet.Letter, alpha alphabet.Alphabet) Seq {
	return Seq{linear.NewSeq(id, b, alpha)}
}

// Clone returns a copy of the Seq.
func (s Seq) Clone() seq.Sequence {
	return Seq{s.Seq.Clone().(*linear.Seq)}
}

// SetDescription sets the Desc of the embedded linear.Seq and parses
// the relevant fields of the description to populate the location, offset
// and strand fields of the sequence annotation.
func (s Seq) SetDescription(d string) error {
	const (
		rangeField  = "range="
		strandField = "strand="
	)

	var (
		start int
		err   error
	)
	for _, f := range strings.Fields(d) {
		switch {
		case strings.HasPrefix(f, rangeField):
			rf := strings.FieldsFunc(f[len(rangeField):], func(r rune) bool { return r == ':' || r == '-' })
			if len(rf) < 1 {
				continue
			}
			s.Loc = Chr(rf[0])
			if len(rf) < 2 {
				continue
			}
			start, err = strconv.Atoi(rf[1])
		case strings.HasPrefix(f, strandField):
			st := f[len(strandField):]
			if len(st) == 0 {
				continue
			}
			if st[0] == '+' {
				s.Strand = seq.Plus
			} else if st[0] == '-' {
				s.Strand = seq.Minus
			}
		}
	}
	s.Offset = feat.OneToZero(start)
	s.Desc = d
	return err
}
