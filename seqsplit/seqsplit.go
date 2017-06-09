// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// seqsplit splits contig sequences that are above a
// minimum cut-off length to generate fragments such that
// each fragment falls in the size range:
//  window ≤ fragment < (2*window).
// This is achieved as follows, calculate:
//  remainder = length % window and
//  quotient = length / window.
// Slice contig from the start till (window * (quotient-1))
// position into (quotient-1) fragments each of size
// window. Slice (last window + remainder) sized fragment
// from the contig starting from position
// (window * (quotient-1)) till the end of the
// contig. Write all fragments to output.
//
// Example: Given a window size of 5kb and a contig of
// size 27582bp, calculate:
//  remainder = 27582 % 5000 = 2582 and
//  quotient = 27582 / 5000 = 5.
// Slice contig from the start till (5000 * (5-1)) position
// into (5-1) fragments each of size 5kb. Get the last
// window+remainder (5000+2582) fragment starting from
// position 20000 till the end of the contig (27582).
package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/seq/sequtils"
)

type fe struct {
	s, e   int
	orient feat.Orientation
	feat.Feature
}

func (f fe) Start() int                    { return f.s }
func (f fe) End() int                      { return f.e }
func (f fe) Len() int                      { return f.e - f.s }
func (f fe) Orientation() feat.Orientation { return f.orient }

type fs []feat.Feature

func (f fs) Features() []feat.Feature { return []feat.Feature(f) }

var (
	inf    = flag.String("in", "", "input contig file name to be fragmented. Defaults to stdin.")
	outf   = flag.String("out", "", "output file name. Defaults to stdout.")
	min    = flag.Int("min", 2500, "minimum sequence length cut-off (bp)")
	window = flag.Int("window", 5000, "sequence window length (bp)")
	help   = flag.Bool("help", false, "help prints this message.")
)

func main() {
	var (
		in, out *os.File
		r       *fasta.Reader
		w       *fasta.Writer
		err     error
	)

	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	t := linear.NewSeq("", nil, alphabet.DNA)
	if *inf == "" {
		r = fasta.NewReader(os.Stdin, t)
	} else if in, err = os.Open(*inf); err != nil {
		log.Fatalf("failed to open %q: %v", *inf, err)
	} else {
		defer in.Close()
		r = fasta.NewReader(in, t)
	}

	if *outf == "" {
		out = os.Stdout
	} else if out, err = os.Create(*outf); err != nil {
		log.Fatalf("failed to create %q: %v", *outf, err)
	}
	defer out.Close()

	w = fasta.NewWriter(out, 60)
	sc := seqio.NewScanner(r)
	for sc.Next() {
		next := sc.Seq().(*linear.Seq)
		curr := linear.NewSeq("", nil, alphabet.DNA)
		startPos, endPos := 0, 0
		switch {
		case len(next.Seq) < *min:
			// Discard contigs below the cut-off size limit.
			continue
		case len(next.Seq) >= 2*(*window):
			remainder := len(next.Seq) % (*window)
			quotient := len(next.Seq) / (*window)
			for i := 0; i < (quotient - 1); i++ {
				startPos = i * (*window)
				endPos = startPos + (*window)
				ff := fs{fe{s: startPos, e: endPos}}
				err := sequtils.Stitch(curr, next, ff)
				if err != nil {
					panic(err)
				}
				// The fragment sequences require new, unique FASTA
				// sequence identifiers. Append the start and end positions
				// of contig sequence to old identifiers and use them as
				// FASTA headers for the fragments.
				curr.Desc = fmt.Sprintf("%v_%v-%v", next.Desc, startPos, endPos)
				if _, err = w.Write(curr); err != nil {
					fmt.Fprintf(os.Stderr, "failed to write window-sized fragment: %v", err)
				}
			}
			ff := fs{fe{s: endPos, e: endPos + (*window) + remainder}}
			err := sequtils.Stitch(curr, next, ff)
			if err != nil {
				panic(err)
			}
			curr.Desc = fmt.Sprintf("%v_%v-%v", next.Desc, endPos, endPos+(*window)+remainder)
			if _, err = w.Write(curr); err != nil {
				fmt.Fprintf(os.Stderr, "failed to write remainder fragment: %v", err)
			}
		default:
			// Contig is of desired size range.
			if _, err = w.Write(next); err != nil {
				fmt.Fprintf(os.Stderr, "failed to write contig: %v", err)
			}
		}
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("failed during read: %v", err)
	}
}
