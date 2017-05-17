// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Split genomic contigs of length above a minimum length into fragments
// of a given window size. If the last remaining fragment is below window
// size join it to the previous fragment, based on:
// https://github.com/tetramerFreqs/Binning/blob/master/tetramer_freqs_esom.pl
package main

import (
	"flag"
	"fmt"
	//	"math"
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
	inf    = flag.String("inf", "test.fna", "input filename")
	outf   = flag.String("outf", "split_test.fna", "output filename")
	min    = flag.Int("min", 2500, "minimum sequence length cut-off (bp)")
	window = flag.Int("window", 5000, "sequence window length (bp)")
	help   = flag.Bool("help", false, "help prints this message.")
)

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	in, err := os.Open(*inf)
	if err != nil {
		fmt.Fprintf(os.Stderr, "open input FASTA file: %v.", err)
		os.Exit(1)
	}
	defer in.Close()
	r := fasta.NewReader(in, linear.NewSeq("", nil, alphabet.DNA))
	out, err := os.Create(*outf)
	if err != nil {
		fmt.Fprintf(os.Stderr, "open output FASTA file: %v.", err)
		os.Exit(1)
	}
	defer out.Close()
	w := fasta.NewWriter(out, 60)
	sc := seqio.NewScanner(r)
	for sc.Next() {
		next := sc.Seq().(*linear.Seq)
 		curr := linear.NewSeq("", nil, alphabet.DNA)
		startPos, endPos := 0, 0
		switch {
		case len(next.Seq) < *min:
			// discard contigs below cut-off size limit
			fmt.Printf("%d bp < %d; discard %s\n", len(next.Seq), *min, next.Name())
		case len(next.Seq) >= 2*(*window):
			//  split contigs and write window-sized fragments
			remainder := len(next.Seq) % *(window)
			quotient := len(next.Seq) / *(window)
			fmt.Printf("l = %d; q = %d; r = %d\n", len(next.Seq), quotient, remainder)
			for i, j := 0, 0;  i < quotient; i, j = i +1, i * (*window) {
				startPos = j
				endPos = startPos + (*window)
				ff := fs{
					fe{s: startPos, e: endPos},
				}
				err := sequtils.Stitch(curr, next, ff)
				if err != nil {
					continue
				}
				// add seq locations to header
				curr.Desc = fmt.Sprintf("%v_%v-%v", next.Desc, startPos, endPos)
				if _, err = w.Write(curr); err != nil {
					fmt.Fprintf(os.Stderr, "failed to write cut fragment :%v", err)
				}
			}
			fmt.Printf("Start: %d, End: %d\n", startPos, endPos)
			// extract and write last remaining fragment
			ff := fs{
				fe{s: endPos, e: endPos + (*window) + remainder},
			}
			err := sequtils.Stitch(curr, next, ff)
				if err != nil {
					continue
				}
			// add seq locations to header
			curr.Desc = fmt.Sprintf("%v_%v-%v", next.Desc, endPos, endPos+(*window)+remainder)
			if _, err = w.Write(curr); err != nil {
				fmt.Fprintf(os.Stderr, "failed to write cut fragment :%v", err)
			}
		default:
			// contig is of desired size range
			if _, err = w.Write(next); err != nil {
				fmt.Fprintf(os.Stderr, "write FASTA record :%v", err)
			}

		}

	}
}
