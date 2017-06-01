// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// seqstats calculates and prints sequence statistics from
// a multi-FASTA DNA sequence file (default stdin). It
// is useful for analyzing metrics of microbial genome
// assemblies or metagenome "bins". It prints: the total
// no. of sequences, assembly size (total length of all
// sequences), Min, Max, Avg and N50.
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

const MaxInt = int(^uint(0) >> 1)

// binStats contains the basename of the file without
// any extension and other reported statistics in bp
// (base pairs)
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
	ctgf = flag.String("ctgf", "", "input contig file, defaults to stdin")
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
	if *ctgf == "" {
		r = fasta.NewReader(os.Stdin, t)
	} else if in, err = os.Open(*ctgf); err != nil {
		log.Fatalf("failed to open %q: %v", *ctgf, err)
		os.Exit(1)
	} else {
		defer in.Close()
		r = fasta.NewReader(in, t)
	}

	sc := seqio.NewScanner(r)
	b.Name = strings.Split(path.Base(*ctgf), ".")[0]
	b.Min = MaxInt

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
	// Print the statistics of the assembly as key-value pairs.
	fmt.Printf("%+v\n", b)
}
