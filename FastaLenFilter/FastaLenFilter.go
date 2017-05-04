// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Filter sequences that are above a length cut-off
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)


var (
	inf     = flag.String("inf", "test.fna", "input filename")
	outf    = flag.String("outf", "test_gt2500bp.fna", "output filename")
	minLen = flag.Int("minLen", 2500, "minimum sequence length cut-off (bp)")
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
			if len(s.Seq) > *minLen {
				_, err := w.Write(s)
				if err != nil {
					fmt.Fprintf(os.Stderr, "write FASTA record: %v", err)
				}
			}
		}
	}
}
