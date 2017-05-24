// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Calculate and print sequence statistics from a multi-FASTA DNA sequence
// file (default stdin), useful for analyzing metrics of microbial genome
// assemblies or metagenome "bins". It prints: the total no. of sequences,
// Assembly size (total length of all sequences), Min, Max, Avg and N50

package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path"
	"sort"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

// binStats contains the basename of the file without any extension and 
// other reported statistics in bp (base pairs)
type binStats struct {
	Name    string // from input filename (empty if stdin)
	totSeqs int
	Size    int
	Min     int
	Max     int
	Avg     int
	N50     int
}

func main() {
	var (
		in      *os.File
		r       *fasta.Reader
		err     error
		b       binStats
		seqlens []int
		inf     = flag.String("inf", "", "input contig file, defaults to stdin")
		help    = flag.Bool("help", false, "help prints this message")
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

	sc := seqio.NewScanner(r)
	b.Name = strings.Split(path.Base(*inf), ".")[0]
	b.Min = 1000000000 // similar: int(^uint(0) >> 1)

	for sc.Next() {
		s := sc.Seq()
		b.totSeqs++
		b.Size += s.Len()
		seqlens = append(seqlens, s.Len())
		if s.Len() < b.Min {
			b.Min = s.Len()
		}
		if s.Len() > b.Max {
			b.Max = s.Len()
		}
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("failed during read: %v", err)
	}
	// sort in descending order of lengths
	sort.Sort(sort.Reverse(sort.IntSlice(seqlens)))
	// csum stores the cumulative sequence lengths
	for csum, i := seqlens[0], 1; i < len(seqlens); i, csum = i+1, seqlens[i]+csum {
		if csum >= (b.Size / 2) {
			b.N50 = seqlens[i]
			break
		}
	}
	b.Avg = b.Size / b.totSeqs
	fmt.Printf("%+v\n", b)
}
