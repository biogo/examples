// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Convert a multiple-sequence alignment in FASTA to
// PHYLIP (sequential) format
package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

func main() {
	var (
		inf            = flag.String("inf", "test.aln", "input FASTA filename")
		outf           = flag.String("outf", "test.phy", "output PHYLIP filename")
		help           = flag.Bool("help", false, "help prints this message.")
		totSeq, alnLen int
	)

	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	// open input (first pass) to get totSeq and alnLen
	in, err := os.Open(*inf)
	if err != nil {
		log.Fatalf("open FASTA file: %v.", err)
		os.Exit(1)
	}
	defer in.Close()
	r := fasta.NewReader(in, linear.NewSeq("", nil, alphabet.Protein))
	out, err := os.Create(*outf)
	if err != nil {
		log.Fatalf("open PHYLIP file: %v.", err)
		os.Exit(1)
	}
	defer out.Close()
	totSeq = 0
	sc := seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq()
		alnLen = s.Len()
		totSeq = totSeq + 1
	}
	alnCounts := fmt.Sprintf("%d %d\n", totSeq, alnLen)
	io.WriteString(out, alnCounts)
	// open input (second pass) to read sequences and write to output
	in, err = os.Open(*inf)
	if err != nil {
		log.Fatalf("open FASTA file: %v.", err)
		os.Exit(1)
	}
	defer in.Close()
	r = fasta.NewReader(in, linear.NewSeq("", nil, alphabet.Protein))
	sc = seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq()
		alnRow := fmt.Sprintf("%s %v\n", s.Name(), s.(*linear.Seq).Seq)
		io.WriteString(out, alnRow)
		if s.Len() != alnLen {
			log.Printf("WARNING: Length of sequence %s is different than %d\n",
				s.Name(), alnLen)
		}
		if len(s.Name()) > 10 {
			log.Printf("WARNING: Sequence ID %s is longer than 10 characters\n",
				s.Name())
		}

	}

}
