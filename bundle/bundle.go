// The bundle program attempts to split multiple FASTA sequence files into
// a collection of multiple FASTA sequence files that are smaller than a
// specified bundle size, omitting sequences below a given size threshold.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

var (
	in     = flag.String("in", "", "Specifies the input filename.")
	cut    = flag.Int("cut", 0, "Specifies the size cut-off for inclusion.")
	bundle = flag.Int("bundle", 20e6, "Specifies the sum of sequence length in a bundle.")
)

func main() {
	flag.Parse()
	if *in == "" {
		flag.Usage()
		os.Exit(1)
	}

	inFile, err := os.Open(*in)
	if err != nil {
		log.Fatalf("failed to open input:%v", err)
	}
	defer inFile.Close()
	*in = filepath.Base(*in)

	sc := seqio.NewScanner(fasta.NewReader(inFile, linear.NewSeq("", nil, alphabet.DNA)))

	var i, size int
	out, err := os.Create(fmt.Sprintf("%s-%d.fa", *in, i))
	for sc.Next() {
		if sc.Seq().Len() < *cut {
			continue
		}
		if size != 0 && size+sc.Seq().Len() > *bundle {
			err = out.Close()
			if err != nil {
				log.Fatalf("failed to close file bundle %d: %v", i, err)
			}
			i++
			size = 0
			out, err = os.Create(fmt.Sprintf("%s-%d.fa", *in, i))
			if err != nil {
				log.Fatalf("failed to open file bundle %d: %v", i, err)
			}
		}
		size += sc.Seq().Len()
		fmt.Fprintf(out, "%60a\n", sc.Seq())
	}
	if sc.Error() != nil {
		log.Fatal(sc.Error())
	}
	err = out.Close()
	if err != nil {
		log.Fatalf("failed to close file bundle %d: %v", i, err)
	}
}
