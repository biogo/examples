// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// fastatophy converts a multiple-sequence alignment in FASTA to
// PHYLIP (sequential) format.

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

var (
	inf  = flag.String("inf", "test.aln", "input FASTA filename")
	outf = flag.String("outf", "test.phy", "output PHYLIP filename")
	help = flag.Bool("help", false, "help prints this message.")
)

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	in, _ := os.Open(*inf)
	defer in.Close()
	r := fasta.NewReader(in, linear.NewSeq("", nil, alphabet.Protein))

	out, _ := os.Create(*outf)
	defer out.Close()

	var n, alnLen int
	n = 0
	sc := seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq()
		alnLen = s.Len()
		n++
	}
	err := sc.Error()
	if err != nil {
		log.Fatalf("failed during read: %v", err)
	}
	fmt.Fprintf(out, "%d %d\n", n, alnLen)

	// Read input FASTA file from start to write sequences to output
	in.Seek(0, io.SeekStart)

	r = fasta.NewReader(in, linear.NewSeq("", nil, alphabet.Protein))
	sc = seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		fmt.Fprintf(out, "%s %v\n", s.Name(), s.Seq)
		if s.Len() != alnLen {
			log.Printf("WARNING: Length of sequence %s is different to %d.\n", s.Name(), alnLen)
		}
		if len(s.Name()) > 10 {
			log.Printf("WARNING: Sequence ID %s is longer than 10 characters.\n", s.Name())
		}

	}

}
