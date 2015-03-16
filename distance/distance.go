// distance reads in a multiple fasta file and compares kmer frequency
// distributions for blocks of each sequece against the sequence average.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/index/kmerindex"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/seq/sequtils"
)

func main() {
	inName := flag.String("in", "", "Filename for input. Defaults to stdin.")
	k := flag.Int("k", 6, "kmer size.")
	chunk := flag.Int("chunk", 1000, "Chunk width.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	var in *fasta.Reader
	if *inName == "" {
		in = fasta.NewReader(os.Stdin, linear.NewSeq("", nil, alphabet.DNA))
	} else if f, err := os.Open(*inName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		in = fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA))
		defer f.Close()
	}

	sub := linear.NewSeq("", nil, alphabet.DNA)
	for {
		s, err := in.Read()
		if err != nil {
			if err != io.EOF {
				fmt.Printf("Error: %v\n", err)
				os.Exit(1)
			}
			return
		}

		ki, err := kmerindex.New(*k, s.(*linear.Seq))
		if err != nil {
			fmt.Printf("Error: %v\n", err)
			os.Exit(1)
		}

		baseLine, ok := ki.NormalisedKmerFrequencies()
		if ok {
			for i := 0; (i+1)**chunk < s.Len(); i++ {
				sequtils.Truncate(sub, s, i**chunk+1, (i+1)**chunk)
				ki, err = kmerindex.New(*k, sub)
				if err != nil {
					fmt.Println(err)
					os.Exit(1)
				}
				if chunkFreqs, ok := ki.NormalisedKmerFrequencies(); ok {
					fmt.Printf("%s\t%d\t%f\n", s.Name(), i**chunk, kmerindex.Distance(baseLine, chunkFreqs))
				}
			}
		}
	}
}
