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
// (base pairs).
type binStats struct {
	name    string // From input filename (empty, if stdin).
	totSeqs int
	size    int
	min     int
	max     int
	avg     float64
	n50     int
}

var (
	ctgf = flag.String("in", "", "input contig file, defaults to stdin")
	help = flag.Bool("help", false, "help prints this message")
)

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	var in *os.File
	var r *fasta.Reader
	var err error
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

	var b binStats
	var seqlens []int
	sc := seqio.NewScanner(r)
	b.name = strings.Split(path.Base(*ctgf), ".")[0]
	b.min = MaxInt

	for sc.Next() {
		s := sc.Seq()
		b.totSeqs++
		b.size += s.Len()
		seqlens = append(seqlens, s.Len())
		if s.Len() < b.min {
			b.min = s.Len()
		}
		if s.Len() > b.max {
			b.max = s.Len()
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
		if csum >= (b.size / 2) {
			b.n50 = seqlens[i]
			break
		}
		csum = seqlens[i] + csum
	}
	b.avg = float64(b.size) / float64(b.totSeqs)
	// Print the statistics of the assembly as key-value pairs.
	fmt.Printf("%+v\n", b)
}
