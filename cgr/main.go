// cgr generates a Chaos Game Representation of a single FASTA sequence.
package main

import (
	"code.google.com/p/biogo/exp/alphabet"
	"code.google.com/p/biogo/exp/seq/linear"
	"code.google.com/p/biogo/exp/seqio/fasta"
	"code.google.com/p/biogo/graphics/color"
	"code.google.com/p/biogo/graphics/kmercolor"
	"code.google.com/p/biogo/index/kmerindex"

	"flag"
	"fmt"
	"image/png"
	"io"
	"os"
)

func main() {
	inName := flag.String("in", "", "Filename for input. Defaults to stdin.")
	outName := flag.String("out", "", "Filename for output. Defaults to stdout.")
	k := flag.Int("k", 6, "kmer size.")
	start := flag.Int("s", 0, "Start site - mandatory parameter > 0.")
	chunk := flag.Int("chunk", 1000, "Chunk width - < 0 indicates sequence to end.")
	desch := flag.Bool("desch", false, "Use diagonal base arrangement described by Deschavanne et al., otherwise use orthogonal arrangement.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	kmerindex.MinKmerLen = *k

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	if *start == 0 {
		fmt.Fprintln(os.Stderr, "Must specify s > 0")
		flag.Usage()
		os.Exit(1)
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

	s, err := in.Read()
	if err != nil {
		if err != io.EOF {
			fmt.Println(err)
			os.Exit(1)
		}
		return
	}
	if *chunk < 0 {
		*chunk = s.Len() - *start - 1
	}

	fmt.Fprintf(os.Stderr, "Indexing %s\n", s.Name())
	ki, err := kmerindex.New(*k, s.(*linear.Seq))
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}

	base := color.HSVA{0, 1, 1, 1}
	cgr := kmercolor.NewCGR(ki, base)
	fmt.Fprintf(os.Stderr, "Painting %s\n", s.Name())
	cgr.Paint(kmercolor.V|kmercolor.H, *desch, *start, *chunk)

	fmt.Fprintf(os.Stderr, "Writing %s\n", s.Name())
	out, err := os.Create(fmt.Sprintf("%s.png", *outName))
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
	}
	png.Encode(out, cgr)
	out.Close()
}
