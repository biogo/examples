// pairwise reads in sequence from two fasta files and reports the
// Euclidian distance between the sequences' kmer distributions.
package main

import (
	"code.google.com/p/biogo/exp/alphabet"
	"code.google.com/p/biogo/exp/seq/linear"
	"code.google.com/p/biogo/exp/seqio/fasta"
	"code.google.com/p/biogo/index/kmerindex"

	"flag"
	"fmt"
	"os"
)

func main() {
	inName1 := flag.String("1", "", "Filename for first input.")
	inName2 := flag.String("2", "", "Filename for second input.")
	k := flag.Int("k", 6, "kmer size.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	var err error

	f1, err := os.Open(*inName1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	}
	defer f1.Close()
	in1 := fasta.NewReader(f1, linear.NewSeq("", nil, alphabet.DNA))

	f2, err := os.Open(*inName2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	}
	defer f2.Close()
	in2 := fasta.NewReader(f2, linear.NewSeq("", nil, alphabet.DNA))

	var (
		kf1, kf2 map[kmerindex.Kmer]float64
		ok       bool
	)

	s1, err := in1.Read()
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	s2, err := in2.Read()
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}

	if ki, err := kmerindex.New(*k, s1.(*linear.Seq)); err != nil {
		fmt.Println(err)
		os.Exit(1)
	} else {
		if kf1, ok = ki.NormalisedKmerFrequencies(); !ok {
			fmt.Printf("Unable to determine Kmer frequences for %s\n", s1.Name())
			os.Exit(1)
		}
	}
	if ki, err := kmerindex.New(*k, s2.(*linear.Seq)); err != nil {
		fmt.Println(err)
		os.Exit(1)
	} else {
		if kf2, ok = ki.NormalisedKmerFrequencies(); !ok {
			fmt.Printf("Unable to determine Kmer frequences for %s\n", s2.Name())
			os.Exit(1)
		}
	}

	fmt.Printf("Kmer distance between %s and %s is %f\n", s1.Name(), s2.Name(), kmerindex.Distance(kf1, kf2))
}
