// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// seqstats calculates and prints sequence statistics from a
// multi-FASTA DNA sequence file (default stdin). It is useful for
// analyzing metrics of microbial genome assemblies or metagenome
// "bins". It prints: the total no. of sequences, assembly size
// (total length of all sequences), Min, Max, Avg and N50
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
	Name    string // from input filename (empty, if stdin)
	totSeqs int
	Size    int
	Min     int
	Max     int
	Avg     int
	N50     int
}

var (
	inf  = flag.String("inf", "", "input contig file, defaults to stdin")
	help = flag.Bool("help", false, "help prints this message")
)

func main() {
	var (
		in      *os.File
		r       *fasta.Reader
		err     error
		b       binStats
		seqlens []int
	)

	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	t := linear.NewSeq("", nil, alphabet.DNA)
	if *inf == "" {
		r = fasta.NewReader(os.Stdin, t)
	} else if in, err = os.Open(*inf); err != nil {
		log.Fatalf("failed to open %q: %v", *inf, err)
		os.Exit(1)
	} else {
		defer in.Close()
		r = fasta.NewReader(in, t)
	}

	sc := seqio.NewScanner(r)
	// Get the basename of the file and remove the extension.
	// Example: "/path/to/infile.fasta" -> "infile.fasta" -> "infile"
	b.Name = strings.Split(path.Base(*inf), ".")[0]
	// assign b.Min to MaxInt64 (1<<63-1 = 9223372036854775807)
	// or MaxInt32 (1<<31-1 = 2147483647)
	b.Min = 1<<63 - 1

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
	// Sort in descending order of sequence length.
	sort.Sort(sort.Reverse(sort.IntSlice(seqlens)))
	// csum stores the cumulative sequence length.
	for i, csum := 1, seqlens[0]; i < len(seqlens); i++ {
		if csum >= (b.Size / 2) {
			b.N50 = seqlens[i]
			break
		}
		csum = seqlens[i] + csum
	}
	b.Avg = b.Size / b.totSeqs
	fmt.Printf("%+v\n", b)
}


