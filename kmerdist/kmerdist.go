// kmerdist performs an analysis on the kmer distribution of a set of sequences.
// It returns summary statistics on the frequencies of kmers in the analysed sequences.
package main

import (
	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/index/kmerindex"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq/linear"

	"flag"
	"fmt"
	"math"
	"os"
	"sort"
)

type Rank []int

func (r Rank) Init() {
	sort.Ints(r)
	if len(r) > 0 {
		for i := range r[1:] {
			r[i+1] += r[i]
		}
	}
}

func (r Rank) At(i int) int {
	if i == 0 {
		return r[0]
	}

	return r[i] - r[i-1]
}

func (r Rank) mass(p float64) int {
	return int(float64(r[len(r)-1]) * p)
}

func (r Rank) Percentile(p float64) float64 {
	switch {
	case r.mass(p) <= r[0]:
		return float64(r.At(0))
	case r.mass(p) >= r[len(r)-1]:
		return float64(r.At(len(r) - 1))
	default:
		min, max := 0, len(r)-1
		for {
			switch mid := (min + max) / 2; {
			case min+1 == max:
				t := float64(r.mass(p))
				l, h := float64(r[min]), float64(r[max])
				return (float64(r.At(min))*(h-t) + float64(r.At(max))*(t-l)) / (h - l)
			case r[mid] == r.mass(p):
				return float64(r.At(mid))
			case r.mass(p) > r[mid]:
				min = mid
			case r.mass(p) < r[mid]:
				max = mid
			}
		}
	}

	panic("cannot reach")
}

func main() {
	inName := flag.String("in", "", "Filename for input. Defaults to stdin.")
	k := flag.Int("k", 6, "kmer size.")
	p := flag.Float64("p", 0.95, "Percentile threshold.")
	fill := flag.Bool("fill", false, "Count NA as 0.")
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

	fk := float64(*k)
	fmt.Printf("ID\tn\tMean\tStDev\tnorm(Mean)\tnorm(StDev)\t95%% percentile\n")
	for {
		s, err := in.Read()
		if err != nil {
			os.Exit(1)
		} else {
			ki, err := kmerindex.New(*k, s.(*linear.Seq))
			if err != nil {
				fmt.Println(err)
				os.Exit(1)
			}
			m, ok := ki.KmerFrequencies()
			if ok {
				r := make(Rank, 0, len(m))
				var n, sumOfSquares, mean, oldmean, kmers float64
				for _, c := range m {
					fc := float64(c)
					kmers += fc
					r = append(r, c)

					// The Method of Provisional Means	
					n++
					mean = oldmean + (fc-oldmean)/n
					sumOfSquares += (fc - oldmean) * (fc - mean)
					oldmean = mean
				}
				r.Init()
				if *fill {
					for n < math.Pow(4, fk) {
						n++
						mean = oldmean * (1 - 1/n)
						sumOfSquares += oldmean * mean
						oldmean = mean
					}
				}
				fl := float64(s.Len())
				stdev := math.Sqrt(sumOfSquares / (n - 1))
				fmt.Printf("%s\t%0.f\t%f\t%f\t%f\t%f\t%f\n",
					s.Name(), n, mean, stdev, mean/fl, stdev/fl, r.Percentile(*p)/kmers)
			}
		}
	}
}
