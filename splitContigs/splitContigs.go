// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Split sequence that is above the min cut-off length to generate fragments
// such that each fragment falls in the range: window < fragment < 2*window
// (i.e. 5 to 10 kb).  To achieve this, split each contig into window-sized
// fragments and add the remainder to the last window-sized fragment.
// Example: given a 14kb long contig,  generate fragments of size 5kb and
// 9kb (rather than 5kb, 5kb and 4kb).
// Add sequence position information to old headers to get new ones.

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

func main() {
	var (
		in, out *os.File
		r       *fasta.Reader
		w       *fasta.Writer
		err     error
		inf     = flag.String("inf", "", "input contig file name to be fragmented. Defaults to stdin.")
		outf    = flag.String("outf", "", "output file name. Defaults to stdout")
		min     = flag.Int("min", 2500, "minimum sequence length cut-off (bp)")
		window  = flag.Int("window", 5000, "sequence window length (bp)")
		help    = flag.Bool("help", false, "help prints this message.")
	)

	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	t := linear.NewSeq("", nil, alphabet.DNA)
	if *inf == "" {
		fmt.Fprintln(os.Stderr, "Reading sequences from stdin.")
		r = fasta.NewReader(os.Stdin, t)
	} else if in, err = os.Open(*inf); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		defer in.Close()
		fmt.Fprintf(os.Stderr, "Reading sequence from `%s'.\n", *inf)
		r = fasta.NewReader(in, t)
	}

	if *outf == "" {
		fmt.Fprintln(os.Stderr, "Writing output to stdout.")
		out = os.Stdout
	} else if out, err = os.Create(*outf); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
	} else {
		fmt.Fprintf(os.Stderr, "Writing output to `%s'.\n", *outf)
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
			// discard contigs below cut-off size limit
			continue
		case len(next.Seq) >= 2*(*window):
			// example: from 27582bp contig, get all 5000 bp sized fragments till
			// 20000 position,  leaving the last fragment + remainder (7582)
			remainder := len(next.Seq) % (*window)
			quotient := len(next.Seq) / (*window)
			for i := 0; i < (quotient - 1); i++ {
				startPos = i * (*window)
				endPos = startPos + (*window)
				ff := fs{
					fe{s: startPos, e: endPos},
				}
				err := sequtils.Stitch(curr, next, ff)
				if err != nil {
					continue
				}
				// add seq locations to header
				curr.Desc = fmt.Sprintf("%v_%v-%v", next.Desc, startPos,
					endPos)
				if _, err = w.Write(curr); err != nil {
					fmt.Fprintf(os.Stderr, "failed to write window-sized fragment: %v", err)
				}
			}
			// example: write last 7582 fragment (5000 + 2582)
			ff := fs{
				fe{s: endPos, e: endPos + (*window) + remainder},
			}
			err := sequtils.Stitch(curr, next, ff)
			if err != nil {
				continue
			}
			curr.Desc = fmt.Sprintf("%v_%v-%v", next.Desc, endPos,
				endPos+(*window)+remainder)
			if _, err = w.Write(curr); err != nil {
				fmt.Fprintf(os.Stderr, "failed to write remainder fragment: %v", err)
			}
		default:
			// contig is of desired size range
			if _, err = w.Write(next); err != nil {
				fmt.Fprintf(os.Stderr, "failed to write contig: %v", err)
			}

		}

	}
}
