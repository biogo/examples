package main

import (
	"flag"
	"fmt"
	"image"
	"image/png"
	"os"

	"github.com/gonum/plot/palette"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/index/kmerindex"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/graphics/kmercolor"
)

func main() {
	var (
		in  *fasta.Reader
		out *os.File
	)

	inName := flag.String("in", "", "Filename for input. Defaults to stdin.")
	outName := flag.String("out", "", "Filename for output. Defaults to stdout.")
	k := flag.Int("k", 6, "kmer size.")
	chunk := flag.Int("chunk", 1000, "Chunk width.")
	height := flag.Int("h", 100, "Rainbow height.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	t := linear.NewSeq("", nil, alphabet.DNA)
	if *inName == "" {
		in = fasta.NewReader(os.Stdin, t)
	} else if f, err := os.Open(*inName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		defer f.Close()
		in = fasta.NewReader(f, t)
	}

	count := 0
	for {
		count++
		if s, err := in.Read(); err != nil {
			os.Exit(1)
		} else {
			if index, err := kmerindex.New(*k, s.(*linear.Seq)); err != nil {
				fmt.Println(err)
				os.Exit(1)
			} else {
				base := palette.HSVA{0, 1, 0, 1}
				rainbow := kmercolor.NewKmerRainbow(image.Rect(0, 0, s.Len() / *chunk, *height), index, base)
				for i := 0; (i+1)**chunk < s.Len(); i++ {
					rainbow.Paint(kmercolor.V, i, *chunk, i, i+1)
				}
				if out, err = os.Create(fmt.Sprintf("%s-%d.png", *outName, count)); err != nil {
					fmt.Fprintf(os.Stderr, "Error: %v.", err)
				}
				png.Encode(out, rainbow)
				out.Close()
			}
		}
	}
}
