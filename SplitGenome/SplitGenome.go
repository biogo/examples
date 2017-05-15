// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Split genomic contigs into fragments of a given size window (winLen).
// Example: if a contig is 14kb long, split into 5 and 9 kb instead of 5, 5 
// and 4 kb fragments, see:
// https://www.ncbi.nlm.nih.gov/pubmed/19698104

package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
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
	inf     = flag.String("inf", "test.fna", "input filename")
	outf    = flag.String("outf", "split_test.fna", "output filename")
	minLen = flag.Int("minLen", 2500, "minimum sequence length cut-off (bp)")
	winLen = flag.Int("winLen", 5000, "sequence window length (bp)")
	help    = flag.Bool("help", false, "help prints this message.")
)

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	f, err := os.Open(*inf)
	if err != nil {
		fmt.Fprintf(os.Stderr, "open input FASTA file: %v.", err)
		os.Exit(1)
	}
	defer f.Close()
	dnaf := fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA))
	of, err := os.Create(*outf)
	if err != nil {
		fmt.Fprintf(os.Stderr, "open output FASTA file: %v.", err)
		os.Exit(1)
	}
	d := linear.NewSeq("", nil, alphabet.DNA)
	res := linear.NewSeq("", nil, alphabet.DNA)
	w := fasta.NewWriter(of, 60)
	defer of.Close()
	for {
		if s, err := dnaf.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				fmt.Fprintf(os.Stderr, "read %v: %v", dnaf, err)
			}
		} else {
			s := s.(*linear.Seq)
			switch {
			case len(s.Seq) >= *minLen && len(s.Seq) < 2*(*winLen):
				if _, err = w.Write(s); err != nil {
					fmt.Fprintf(os.Stderr, "write FASTA record :%v", err)
				}
			case len(s.Seq) >= 2*(*winLen):
				for i := 0; i < len(s.Seq); i = i + *winLen {
					j := len(s.Seq) - i
					if j < *winLen {
						prev := d
						ff := fs{
							fe{s: i, e: j + i + *winLen},
						}
						if err := sequtils.Stitch(res, s, ff); err == nil {
							res.Desc = fmt.Sprintf("%v_clean_%v", s.Desc, i)
							sequtils.Join(prev, res, 2)
						}
						prev.Desc = fmt.Sprintf("%v_remainder_%v", s.Desc,
							len(prev.Seq))
						if _, err = w.Write(prev); err != nil {
							fmt.Fprintf(os.Stderr, "write FASTA record :%v", err)
						}
					} else {
						ff := fs{
							fe{s: i, e: i + *winLen},
						}
						if err := sequtils.Stitch(d, s, ff); err == nil {
							d.Desc = fmt.Sprintf("%v_clean_%v", s.Desc, i)
							if _, err = w.Write(d); err != nil {
								fmt.Fprintf(os.Stderr, "write FASTA record :%v", err)
							}
						}
					}
				}
			case len(s.Seq) < *minLen:
				fmt.Printf("%d bp < %d; discard %s\n", len(s.Seq), *minLen, s.Name())
			}

		}

	}
}
