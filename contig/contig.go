// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package contig provides storage and representation of sequences constructed
// from subsequence contigs.
package contig

import (
	"errors"
	"fmt"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/sequtils"
	"github.com/biogo/biogo/util"
	"github.com/biogo/store/step"
)

type seqStep struct {
	seq.Sequence
}

// Equal returns a boolean indicating equality between the receiver
// and the parameter. Equal will panic if the parameter does not satisfy
// the seq.Sequence interface or is an ambig.
func (s seqStep) Equal(e step.Equaler) bool {
	switch e := e.(type) {
	case ambig:
		return false
	case seqStep:
		if s.Sequence == e.Sequence {
			return true
		}
		return false
	}
	panic("contig: non-sequence types not handled")
}

func (s seqStep) String() string { return fmt.Sprint(s.Sequence) }

type ambig alphabet.Letter

func (a ambig) Equal(e step.Equaler) bool   { return a == e.(ambig) }
func (a ambig) Format(fs fmt.State, _ rune) { fs.Write([]byte{byte(a)}) }

type Contig struct {
	*seq.Annotation
	vector *step.Vector
}

// New returns a new super contig sequence spanning the positions [0, l) and
// using the provided alphabet's ambiguous letter as the step ground state.
func New(id string, l int, a alphabet.Alphabet) (*Contig, error) {
	v, err := step.New(0, l, ambig(a.Ambiguous()))
	if err != nil {
		return nil, err
	}
	return &Contig{
		vector:     v,
		Annotation: &seq.Annotation{ID: id, Alpha: a},
	}, nil
}

// Relaxed sets the Contig's length restriction relaxation to the boolean r.
func (c *Contig) Relaxed(r bool) { c.vector.Relaxed = r }

// IsRelaxed returns whether the Contig allows insertion of contigs outside its length.
func (c *Contig) IsRelaxed() bool { return c.vector.Relaxed }

// Joiner returns the ground state of the Contig.
func (c *Contig) Joiner() alphabet.Letter { return alphabet.Letter(c.vector.Zero.(ambig)) }

// Insert adds a sequence to the Contig. The sequence's alphabet must match the Contig's
// alphabet. If the Contig is not relaxed an insertion beyond the range of the contig will
// return an out of range error.
func (c *Contig) Insert(s seq.Sequence) error {
	if s.Alphabet() != c.Alphabet() {
		return errors.New("contig: alphabet mismatch")
	}
	if !c.vector.Relaxed && s.Start() < 0 || s.End() > s.End() {
		return errors.New("contig: sequence out of range")
	}
	c.insert(s)
	return nil
}

func (c *Contig) insert(s seq.Sequence) {
	c.vector.SetRange(s.Start(), s.End(), seqStep{s})
}

// Start returns the Start position of the Contig.
func (c *Contig) Start() int { return c.vector.Start() }

// End returns the End position of the Contig.
func (c *Contig) End() int { return c.vector.End() }

// Len returns the length of the Contig.
func (c *Contig) Len() int { return c.vector.Len() }

// At returns the letter at position i of the Contig. At will panic if i is outside
// the range of the Contig.
func (c *Contig) At(i int) alphabet.QLetter {
	l, err := c.vector.At(i)
	if err != nil {
		panic("err")
	}
	switch l := l.(type) {
	case ambig:
		return alphabet.QLetter{L: alphabet.Letter(l), Q: seq.DefaultQphred}
	case seq.Sequence:
		return l.At(i)
	}
	panic("contig: non-seq type not handled")
}

// Set sets the letter at postion i to l. If no sequence is present at the specified
// position, Set is a no-op on the Contig and returns a non-nil error.
func (c *Contig) Set(i int, l alphabet.QLetter) error {
	vs, err := c.vector.At(i)
	if err != nil {
		return err
	}
	switch vs := vs.(type) {
	case ambig:
		return errors.New("contig: no sequence at specified position")
	case seq.Sequence:
		vs.Set(i, l)
		return nil
	}
	panic("contig: non-seq type not handled")
}

// RevComp reverse complements the Contig and its contained sequences.
func (c *Contig) RevComp() {
	v, _ := step.New(0, c.vector.Len(), c.vector.Zero)
	v.Relaxed = c.vector.Relaxed
	c.vector.Do(func(start, end int, e step.Equaler) {
		if e, ok := e.(seqStep); ok {
			e.RevComp()
			e.SetOffset(c.Len() - e.End())
			v.SetRange(e.Start(), e.End(), e)
		}
	})
	c.vector = v
	c.Strand = -c.Strand
	fmt.Println()
}

// Reverse reverses the Contig and its contained sequences.
func (c *Contig) Reverse() {
	v, _ := step.New(0, c.vector.Len(), c.vector.Zero)
	v.Relaxed = c.vector.Relaxed
	c.vector.Do(func(start, end int, e step.Equaler) {
		if e, ok := e.(seqStep); ok {
			e.Reverse()
			e.SetOffset(c.Len() - e.End())
			v.SetRange(e.Start(), e.End(), e)
		}
	})
	c.vector = v
	c.Strand = seq.None
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// Format is a fmt.Formatter helper. It provides support for the %v (with go syntax
// representation), %s and %a (FASTA output). Note that %v representation takes no
// account of overlapping inserted sequences and so will be misleading if overlaps
// exist.
func (c *Contig) Format(fs fmt.State, cr rune) {
	if c == nil {
		fmt.Fprint(fs, "<nil>")
		return
	}
	var (
		w, _   = fs.Width()
		p, pOk = fs.Precision()
		limit  = -1
	)
	if pOk {
		limit = min(p, c.Len())
	} else {
		p = c.Len()
	}
	switch cr {
	case 'v':
		if fs.Flag('#') {
			fmt.Fprintf(fs, "&%#v", *c)
			return
		} else {
			fmt.Fprintf(fs, "%s", c.vector)
		}
		return
	case 's':
		if !fs.Flag('-') {
			fmt.Fprintf(fs, "%q ", c.ID)
		}
		w = 0
	case 'a':
		fmt.Fprintf(fs, "%c%s", '>', c.ID)
		if c.Desc != "" {
			fmt.Fprintf(fs, " %s", c.Desc)
		}
		fmt.Fprintln(fs)
	default:
		fmt.Fprintf(fs, "%%!%c(contig.Contig=%.10s)", cr, c)
		return
	}
	lw := util.NewWrapper(fs, w, limit)
	c.vector.DoRange(0, p,
		func(start, end int, e step.Equaler) {
			switch e := e.(type) {
			case seqStep:
				if e.Start() != start || e.End() != end {
					se := e.New()
					sequtils.Truncate(se, e, start, end)
					fmt.Fprintf(lw, "%-s", se)
					break
				}
				fmt.Fprintf(lw, "%-s", e.Sequence)
			case ambig:
				eb := []byte{byte(e)}
				for i := start; i < end; i++ {
					lw.Write(eb)
				}
			}
		},
	)
	if pOk && p < c.Len() {
		fmt.Fprint(fs, "...")
	}
}
