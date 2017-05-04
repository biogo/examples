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
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)
var (
	inf     = flag.String("inf", "test.aln", "input FASTA filename")
	outf    = flag.String("outf", "test.phy", "output PHYLIP filename")
	help    = flag.Bool("help", false, "help prints this message.")
)

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	f, err := os.Open("test.aln")
	var (
		alnLen int // alignment length
		totSeq int // Total number of sequences
	)
	if err != nil {
		fmt.Fprintf(os.Stderr, "open FASTA file: %v.", err)
		os.Exit(1)
	}
	defer f.Close()
	dnaf := fasta.NewReader(f, linear.NewSeq("", nil, alphabet.Protein))
	of, err := os.Create("test.phy")
	if err != nil {
		fmt.Fprintf(os.Stderr, "open PHYLIP file: %v.", err)
		os.Exit(1)
	}
	defer of.Close()
	totSeq = 0
	for {
		if s, err := dnaf.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				fmt.Printf("Failed to read %q: %s", dnaf, err)
			}
		} else {
			s := s.(*linear.Seq)
			alnLen = len(s.Seq)
			totSeq = totSeq + 1
		}
	}
	f, err = os.Open(*inf)
	if err != nil {
		fmt.Fprintf(os.Stderr, "open FASTA file: %v.", err)
		os.Exit(1)
	}
	defer f.Close()
	dnaf = fasta.NewReader(f, linear.NewSeq("", nil, alphabet.Protein))
	alnCounts := fmt.Sprintf("%d %d\n", totSeq, alnLen)
	io.WriteString(of, alnCounts)
	for {
		if s, err := dnaf.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				fmt.Printf("Failed to read %q: %s", dnaf, err)
			}
		} else {
			s := s.(*linear.Seq)
			if len(s.Seq) != alnLen {
				log.Printf("WARNING: Length of sequence %s is different than %d\n",
					s.Name(), alnLen)
			}
			if len(s.Name()) > 10 {
				log.Printf("WARNING: Sequence ID %s is longer than 10 characters\n",
					s.Name())
			}
			alnRow := fmt.Sprintf("%s %s\n", s.Name(), s.Seq)
			io.WriteString(of, alnRow)

		}
	}

}
