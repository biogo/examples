// Copyright ©2017 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// seqlen filters sequences that are above a length cut-off.
package main

import (
	"flag"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

// Flags assign input, output file names and minimum length cut-off to their default
// values and also provide a way to change them easily  in the command-line.
// An example using flags for file names, with default sequence length cut-off of 2500 bp.
// seqlen -inf infile.fasta -outf outfile.fasta
// An example using stdin, stdout, with sequence length cut-off of 200 bp.
// seqlen < infile.fasta > outfile.fasta -min 200
var (
	inf  = flag.String("inf", "", "input contig file name to be fragmented. Defaults to stdin.")
	outf = flag.String("outf", "", "output file name. Defaults to stdout")
	min  = flag.Int("min", 2500, "minimum sequence length cut-off (bp)")
	help = flag.Bool("help", false, "help prints this message.")
)

func main() {
	var (
		in, out *os.File
		r       *fasta.Reader
		w       *fasta.Writer
		err     error
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
	} else {
		defer in.Close()
		r = fasta.NewReader(in, t)
	}

	if *outf == "" {
		out = os.Stdout
	} else if out, err = os.Create(*outf); err != nil {
		log.Fatalf("failed to open %q: %v", *outf, err)
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

