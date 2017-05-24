// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Filter sequences that are above a length cut-off
package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

func main() {
	var (
		in, out *os.File
		r       *fasta.Reader
		w       *fasta.Writer
		err     error
		inf     = flag.String("inf", "", "input contig file name to be fragmented. Defaults to stdin.")
		outf    = flag.String("outf", "", "output file name. Defaults to stdout")
		min     = flag.Int("minLen", 2500, "minimum sequence length cut-off (bp)")
		help    = flag.Bool("help", false, "help prints this message.")
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

	if *outf == "" {
		fmt.Fprintln(os.Stderr, "Writing output to stdout.")
		out = os.Stdout
	} else if out, err = os.Create(*outf); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
	} else {
		fmt.Fprintf(os.Stderr, "Writing output to `%s'.\n", *outf)
	}
	defer out.Close()
	w = fasta.NewWriter(out, 60)
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
