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
	curr := linear.NewSeq("", nil, alphabet.DNA)
	prev := linear.NewSeq("", nil, alphabet.DNA)
	w := fasta.NewWriter(out, 60)
	sc := seqio.NewScanner(r)
	for sc.Next() {
		next := sc.Seq().(*linear.Seq)
		switch {
		default:
			// contig is of the desired size range, write without any further modification
			if _, err = w.Write(next); err != nil {
				fmt.Fprintf(os.Stderr, "write FASTA record :%v", err)
			}
		case len(next.Seq) >= 2*(*window):
			// contig is over twice the window size, split it into window-sized fragments
			for i := 0; i < len(next.Seq); i = i + *window {
				j := len(next.Seq) - i
				if j > *window {
					// remaining fragment size is above window size
					ff := fs{
						fe{s: i, e: i + *window},
					}
					if err := sequtils.Stitch(curr, next, ff); err == nil {
						curr.Desc = fmt.Sprintf("%v_clean_%v", next.Desc, i)
						if _, err = w.Write(curr); err != nil {
							fmt.Fprintf(os.Stderr, "write FASTA record :%v", err)
						}
					}
				} else {
					// remaining fragment size is below window size,
					// join it to the previous fragment
					ff := fs{
						fe{s: i, e: j + i + *window},
					}
					if err := sequtils.Stitch(prev, next, ff); err == nil {
						sequtils.Join(curr, prev, 2)
					}
					curr.Desc = fmt.Sprintf("%v_remainder_%v", next.Desc, len(curr.Seq))
					if _, err = w.Write(curr); err != nil {
						fmt.Fprintf(os.Stderr, "write FASTA record :%v", err)
					}
				}
			}
		case len(next.Seq) < *min:
			// discard contigs if they are below the minimum size cutoff
			fmt.Printf("%d bp < %d; discard %s\n", len(next.Seq), *min, next.Name())
		}

	}
}
