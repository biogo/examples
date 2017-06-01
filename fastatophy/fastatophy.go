// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// fastatophy converts a multiple-sequence alignment in
// FASTA to PHYLIP (sequential) format.
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

// padding returns requested number of spaces as a
// concatenated string.
func padding(n int) (spaces string) {
	space := " "
	for i := 0; i < n; i++ {
		spaces += space
	}
	return spaces
}

var (
	inf  = flag.String("inf", "test.aln", "input FASTA filename")
	outf = flag.String("outf", "test.phy", "output PHYLIP filename")
	help = flag.Bool("help", false, "help prints this message")
)

func main() {
	var (
		n, alnLen  int
		strictName string
	)

	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	in, err := os.Open(*inf)
	if err != nil {
		log.Fatalf("failed to open FASTA file: %v", err)
	}
	defer in.Close()
	r := fasta.NewReader(in, linear.NewSeq("", nil, alphabet.Protein))

	out, err := os.Create(*outf)
	if err != nil {
		log.Fatalf("failed to open PHYLIP file: %v", err)
	}
	defer out.Close()

	// Read all FASTA records to get total number of sequences
	// (n) and length of each sequence. alnLen stores the
	// sequence length of the last record.
	n = 0
	sc := seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq()
		alnLen = s.Len()
		n++
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("failed during read: %v", err)
	}
	// Write the header section consisting of dimensions of
	// the alignment to the PHYLIP file.
	fmt.Fprintf(out, "%d %d\n", n, alnLen)

	// Reinitialize to read from the start of the FASTA file
	// and write the alignment section to the PHYLIP file.
	if _, err := in.Seek(0, io.SeekStart); err != nil {
		log.Fatalf("seek failed: %v", err)
	}
	r = fasta.NewReader(in, linear.NewSeq("", nil, alphabet.Protein))
	sc = seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		// Assert that each sequence in the multiple-sequence
		// alignment is of equal length.
		if s.Len() != alnLen {
			log.Printf("Identifier: %s length differs from %d\n", s.Name(), alnLen)
		}
		// Sequence identifiers must be exactly 10 characters in
		// "strict" PHYLIP format, truncate to first 10 characters
		// if identifiers are longer, otherwise pad them with
		// spaces.
		if len(s.Name()) > 10 {
			strictName = s.Name()[:10]
			log.Printf("Identifier: %s was truncated to 10 characters\n", s.Name())
		} else {
			strictName = s.Name() + padding(10-len(s.Name()))
		}
		fmt.Fprintf(out, "%s %v\n", strictName, s.Seq)
	}
}
