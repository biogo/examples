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

var (
	inf  = flag.String("in", "", "input FASTA filename")
	outf = flag.String("out", "", "output PHYLIP filename")
	help = flag.Bool("help", false, "help prints this message")
)

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	t := linear.NewSeq("", nil, alphabet.Protein)
	if *inf == "" {
		flag.Usage()
		os.Exit(1)
	} 
	
	var in *os.File
	var r *fasta.Reader
	var err error
	in, err = os.Open(*inf)
	if err != nil {
		log.Fatalf("failed to open FASTA file: %v", err)
	}
	defer in.Close()
	r = fasta.NewReader(in, t)

	var out *os.File
	if *outf == "" {
		flag.Usage()
	} else if out, err = os.Create(*outf); err != nil {
		log.Fatalf("failed to open PHYLIP file: %v", err)
	} else {
		defer out.Close()
	}

	// Read all FASTA records to get total number of sequences
	// (n) and length of each sequence (seqlens).
	var n int
	var seqlens []int

	sc := seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq()
		seqlens = append(seqlens, s.Len())
		// Assert that each sequence in the multiple-sequence
		// alignment is of equal length.
		if n > 0 {
			if s.Len() != seqlens[n-1] {
				log.Printf("%s length (%d) differs from previous sequence (%d) \n", s.Name(), s.Len(), seqlens[n-1])
			}
		}
		n++
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("failed during first read: %v", err)
	}

	// Write the header section consisting of dimensions of
	// the alignment to the PHYLIP file.
	fmt.Fprintf(out, "%d %d\n", n, seqlens[n-1])

	// Reinitialize to read from the start of the FASTA file
	// and write the alignment section to the PHYLIP file.
	_, err = in.Seek(0, io.SeekStart)
	if err != nil {
		log.Fatalf("seek failed: %v", err)
	}
	r = fasta.NewReader(in, t)
	sc = seqio.NewScanner(r)
	var strictName string
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		// Sequence identifiers must be exactly 10 characters in
		// "strict" PHYLIP format, truncate to first 10 characters
		// if identifiers are longer, otherwise pad them with
		// spaces.
		if len(s.Name()) > 10 {
			strictName = s.Name()[:10]
			log.Printf("Identifier: %s was truncated to 10 characters\n", s.Name())
		} else {
			const padding = "          " // Ten spaces.
			strictName = s.Name() + padding[:10-len(s.Name())]
		}
		fmt.Fprintf(out, "%s %v\n", strictName, s.Seq)
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("failed during second read: %v", err)
	}
}
