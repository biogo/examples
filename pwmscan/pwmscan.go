// pwmscan performs a position weight matrix scan of a set of sequences to
// search for a motif.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/pwm"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/seq/multi"
)

func main() {
	var (
		in, min *fasta.Reader
		mf      *os.File
		matin   *bufio.Reader
		align   *multi.Multi
		out     *gff.Writer
		err     error
	)

	inName := flag.String("in", "", "Filename for input. Defaults to stdin.")
	matName := flag.String("mat", "", "Filename for matrix/alignment input.")
	num := flag.Bool("num", false, "Use numerical description rather than sequence.")
	outName := flag.String("out", "", "Filename for output. Defaults to stdout.")
	precision := flag.Int("prec", 6, "Precision for floating point output.")
	minScore := flag.Float64("score", 0.9, "Minimum score for a hit.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	if *matName == "" {
		flag.Usage()
		os.Exit(1)
	}

	matrix := [][]float64{}
	if *num {
		if mf, err = os.Open(*matName); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
			os.Exit(1)
		} else {
			matin = bufio.NewReader(mf)
		}
		defer mf.Close()

		for {
			line, err := matin.ReadBytes('\n')
			if err != nil {
				break
			}
			if line[len(line)-1] == '\n' {
				line = line[:len(line)-1]
			}
			fields := strings.Split(string(line), "\t")
			if len(fields) < 4 {
				break
			}
			matrix = append(matrix, make([]float64, 0, 4))
			for _, s := range fields {
				if f, err := strconv.ParseFloat(s, 64); err != nil {
					fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
					os.Exit(1)
				} else {
					matrix[len(matrix)-1] = append(matrix[len(matrix)-1], f)
				}
			}
		}
	} else {
		mr, err := os.Open(*matName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
			os.Exit(1)
		}
		min = fasta.NewReader(mr, linear.NewSeq("", nil, alphabet.DNA))
		for {
			s, err := min.Read()
			if err != nil {
				if err != io.EOF {
					fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
					os.Exit(1)
				}
				break
			}
			align.Add(s)
		}
		mr.Close()

		for i := 0; i < align.Len(); i++ {
			matrix[i] = make([]float64, 4)
			for _, v := range align.Column(i, true) {
				if base := alphabet.DNA.IndexOf(v); base >= 0 {
					matrix[i][base]++
				}
			}
		}
	}
	wm := pwm.New(matrix)
	wm.Format = fmt.Sprintf("%%.%de", *precision)

	var r io.ReadCloser
	if *inName == "" {
		r = os.Stdin
	} else if r, err = os.Open(*inName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
		os.Exit(1)
	} else {
		defer r.Close()
	}
	in = fasta.NewReader(r, linear.NewSeq("", nil, alphabet.DNA))

	if *matName == "" {
		flag.Usage()
		os.Exit(1)
	}

	if *outName == "" {
		out = gff.NewWriter(os.Stdout, 60, true)
	} else if f, err := os.Create(*outName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
	} else {
		defer f.Close()
		buf := bufio.NewWriter(f)
		defer buf.Flush()
		out = gff.NewWriter(buf, 60, true)
	}
	out.Precision = 2

	for {
		if s, err := in.Read(); err != nil {
			break
		} else {
			fmt.Fprintf(os.Stderr, "Working on: %s %s\n", s.Name(), s.Description())

			res := wm.Search(s.(*linear.Seq), s.Start(), s.End(), *minScore)
			if len(res) == 1 {
				fmt.Fprintf(os.Stderr, "... found %d match.\n", len(res))
			} else {
				fmt.Fprintf(os.Stderr, "... found %d matches.\n", len(res))
			}
			if len(res) > 0 {
				out.WriteMetaData(gff.Sequence{s.Name(), s.Alphabet().Moltype()})
			}
			for _, r := range res {
				m := r.(*pwm.Feature)
				out.Write(&gff.Feature{
					SeqName:    s.Name(),
					Source:     "pwmscan",
					Feature:    "match",
					FeatStart:  m.MotifStart,
					FeatEnd:    m.MotifEnd,
					FeatScore:  &m.MotifScore,
					FeatStrand: seq.Strand(m.MotifOrient),
					FeatFrame:  gff.NoFrame,
					FeatAttributes: gff.Attributes{
						gff.Attribute{
							Tag:   "Motif",
							Value: fmt.Sprintf("%-v", m.MotifSeq),
						},
						gff.Attribute{
							Tag:   "p",
							Value: fmt.Sprintf("%.*f", *precision, m.MotifProb),
						},
					},
				})
			}
		}
	}
}
