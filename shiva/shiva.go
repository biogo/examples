package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"runtime/pprof"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

func main() {
	var (
		in      *fasta.Reader
		out     *fasta.Writer
		err     error
		profile *os.File
	)

	inName := flag.String("in", "", "Filename for input. Defaults to stdin.")
	outName := flag.String("out", "", "Filename for output. Defaults to stdout.")
	size := flag.Int("size", 40, "Fragment size.")
	width := flag.Int("width", 60, "Fasta output width.")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to this file.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	if *cpuprofile != "" {
		if profile, err = os.Create(*cpuprofile); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.", err)
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Writing CPU profile data to %s\n", *cpuprofile)
		pprof.StartCPUProfile(profile)
		defer pprof.StopCPUProfile()
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

	if *outName == "" {
		out = fasta.NewWriter(os.Stdout, *width)
	} else if f, err := os.Create(*outName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		defer f.Close()
		buf := bufio.NewWriter(f)
		defer buf.Flush()
		out = fasta.NewWriter(buf, *width)
	}

	var trunc *linear.Seq
	for {
		s, err := in.Read()
		if err != nil {
			break
		}
		length := s.Len()
		li := s.(*linear.Seq)
		trunc.ID = li.ID
		switch {
		case length >= 20 && length <= 85:
			t.Seq = li.Seq[5:]
			out.Write(t)
		case length > 85:
			for start := 0; start+*size <= length; start += *size {
				t.Seq = li.Seq[start : start+*size]
				out.Write(t)
			}
		}
	}
}
