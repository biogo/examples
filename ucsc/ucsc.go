// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// ucsc is an example client of biogo.examples/ucsc/ucsc. I reads a set of
// UCSC exported FASTA sequences and parses the location information in the
// sequence description into the sequence metadata, converting the 1-based
// UCSC position information into 0-based half-open used by bíogo. It then
// prints out the FASTA, preceded by a summary of the location.
package main

import (
	"fmt"
	"io"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"

	"github.com/biogo/examples/ucsc/ucsc"
)

var fa = `>hg19_dna range=chr18:78016000-78016181 5'pad=0 3'pad=0 strand=+ repeatMasking=none
AGAGGGAGGATTATTATAATATTGGAATAAAGAGTAATTGCTATCAACTA
ATGATTAATGATATTCATATATAATCATGTCTAAGATCTATATCTGGTAT
AACTATTCTTGTTTTATATTTTATTGGAGTGGAACAGCTCATGTCCTCGG
TCTCTTGCCTCGGCAAAGATTAGATTAGGGTT
>hg19_dna range=chr18:78015995-78016181 5'pad=5 3'pad=0 strand=+ repeatMasking=none
ATTATAGAGGGAGGATTATTATAATATTGGAATAAAGAGTAATTGCTATC
AACTAATGATTAATGATATTCATATATAATCATGTCTAAGATCTATATCT
GGTATAACTATTCTTGTTTTATATTTTATTGGAGTGGAACAGCTCATGTC
CTCGGTCTCTTGCCTCGGCAAAGATTAGATTAGGGTT
>hg19_dna range=chr18:78016000-78016181 5'pad=0 3'pad=0 strand=- repeatMasking=none
AACCCTAATCTAATCTTTGCCGAGGCAAGAGACCGAGGACATGAGCTGTT
CCACTCCAATAAAATATAAAACAAGAATAGTTATACCAGATATAGATCTT
AGACATGATTATATATGAATATCATTAATCATTAGTTGATAGCAATTACT
CTTTATTCCAATATTATAATAATCCTCCCTCT
>hg19_dna range=chr18:78016000-78016186 5'pad=5 3'pad=0 strand=- repeatMasking=none
CACCTAACCCTAATCTAATCTTTGCCGAGGCAAGAGACCGAGGACATGAG
CTGTTCCACTCCAATAAAATATAAAACAAGAATAGTTATACCAGATATAG
ATCTTAGACATGATTATATATGAATATCATTAATCATTAGTTGATAGCAA
TTACTCTTTATTCCAATATTATAATAATCCTCCCTCT
`

func main() {
	buf := strings.NewReader(fa)
	r := fasta.NewReader(buf, ucsc.NewSeq("", nil, alphabet.DNA))
	for {
		s, err := r.Read()
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		fmt.Printf("\nChr:%s Start:%d End:%d Len:%d Strand:%v\n\n%60a\n",
			s.Location(), s.Start(), s.End(), s.Len(), seq.Strand(s.(feat.Orienter).Orientation()), s)
	}
}
