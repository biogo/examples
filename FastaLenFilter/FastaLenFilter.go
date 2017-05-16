// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Filter sequences that are above a length cut-off
package main

import (
	"flag"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

var (
	inf  = flag.String("inf", "test.fna", "input filename")
	outf = flag.String("outf", "test_gt2500bp.fna", "output filename")
	min  = flag.Int("minLen", 2500, "minimum sequence length cut-off (bp)")
	help = flag.Bool("help", false, "help prints this message.")
)

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	in, err := os.Open(*inf)
	if err != nil {
		log.Fatalf("open input FASTA file: %v.", err)
		os.Exit(1)
	}
	defer in.Close()
	r := fasta.NewReader(in, linear.NewSeq("", nil, alphabet.DNA))
	out, err := os.Create(*outf)
	if err != nil {
		log.Fatalf("open output FASTA file: %v.", err)
		os.Exit(1)
	}
	w := fasta.NewWriter(out, 60)
	defer out.Close()
	sc := seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq()
		if s.Len() > *min {
			_, err := w.Write(s)
			if err != nil {
				log.Fatalf("failed to write sequence %q: %v", s.Name(), err)
			}
		}
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("failed during read: %v", err)
	}
}
