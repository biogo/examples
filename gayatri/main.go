package main

import (
	"code.google.com/p/biogo.matrix"
	"code.google.com/p/biogo/exp/seqio/fasta"
	"code.google.com/p/biogo/index/kmerindex"

	"code.google.com/p/biogo/exp/alphabet"
	"code.google.com/p/biogo/exp/seq/linear"
	"flag"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"
	"sort"
	"time"
)

func main() {
	var (
		in           *fasta.Reader
		out, profile *os.File
		err          error
	)

	inName := flag.String("in", "", "Filename for input to be factorised. Defaults to stdin.")
	outName := flag.String("out", "", "Filename for output. Defaults to stdout.")
	k := flag.Int("k", 8, "kmer size to use.")
	cat := flag.Int("cat", 5, "number of categories.")
	iter := flag.Int("i", 1000, "iterations.")
	limit := flag.Int("time", 10, "time limit for NMF.")
	lo := flag.Int("lo", 1, "minimum number of kmer frequency to use in NMF.")
	hi := flag.Float64("hi", 0.5, "maximum proportion of kmer representation to use in NMF.")
	tol := flag.Float64("tol", 0.001, "tolerance for NMF.")
	seed := flag.Int64("seed", -1, "seed for random number generator (-1 uses system clock).")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to this file.")
	help := flag.Bool("help", false, "print this usage message.")

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
		fmt.Fprintln(os.Stderr, "Reading sequences from stdin.")
		in = fasta.NewReader(os.Stdin, t)
	} else if f, err := os.Open(*inName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		defer f.Close()
		fmt.Fprintf(os.Stderr, "Reading sequence from `%s'.\n", *inName)
		in = fasta.NewReader(f, t)
	}

	if *outName == "" {
		fmt.Fprintln(os.Stderr, "Writing output to stdout.")
		out = os.Stdout
	} else if out, err = os.Create(*outName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
	} else {
		fmt.Fprintf(os.Stderr, "Writing output to `%s'.\n", *outName)
	}
	defer out.Close()

	totalkmers := make(map[kmerindex.Kmer]float64)
	kmerlists := make([]map[kmerindex.Kmer]float64, 0)
	seqTable := make([]string, 0)

	for {
		if s, err := in.Read(); err != nil {
			break
		} else {
			var freqs map[kmerindex.Kmer]float64
			if kindex, err := kmerindex.New(*k, s.(*linear.Seq)); err != nil {
				fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
				os.Exit(1)
			} else {
				freqs, _ = kindex.NormalisedKmerFrequencies()
				kmerlists = append(kmerlists, freqs)
				for kmer, freq := range freqs {
					totalkmers[kmer] += freq
				}
			}
			seqTable = append(seqTable, string(s.Name()))
		}
	}

	kmerArray := make([][]float64, 0)
	kmerTable := make([]kmerindex.Kmer, 0)

	for kmer, _ := range totalkmers {
		var count int
		for _, kmerlist := range kmerlists {
			if kmerlist[kmer] > 0 {
				count++
			}
		}
		if count < *lo || float64(count)/float64(len(kmerlists)) > *hi {
			continue
		}
		row := make([]float64, len(kmerlists))
		for i, kmerlist := range kmerlists {
			row[i] = float64(kmerlist[kmer])
		}
		kmerArray = append(kmerArray, row)
		kmerTable = append(kmerTable, kmer)
	}

	kMat, err := matrix.NewDense(kmerArray)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	}

	f := func(i, j int, v float64) float64 {
		if kMat.At(i, j) != 0 {
			return 1
		}
		return 0
	}
	nonZero := kMat.Apply(f, kMat).Sum()

	r, c := kMat.Dims()
	density := nonZero / float64(r*c)

	if *seed == -1 {
		*seed = time.Now().UnixNano()
	}
	fmt.Fprintf(os.Stderr, "Using %v as random seed.\n", *seed)
	rand.Seed(*seed)

	rows, cols := kMat.Dims()

	posNorm := func() float64 { return math.Abs(rand.NormFloat64()) }

	Wo, err := matrix.FuncDense(rows, *cat, 1, posNorm)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	}

	Ho, err := matrix.FuncDense(*cat, cols, 1, posNorm)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	}

	fmt.Fprintf(os.Stderr, "Dimensions of Kmer matrix = (%v, %v)\nDensity = %.3f %%\n%v\n", r, c, (density)*100, kMat)

	W, H, ok := matrix.Factors(kMat, Wo, Ho, *tol, *iter, time.Duration(*limit)*1e9)

	fmt.Fprintf(os.Stderr, "norm(H) = %v norm(W) = %v\n\nFinished = %v\n\n", H.Norm(matrix.Fro), W.Norm(matrix.Fro), ok)

	printFeature(out, kMat, W, H, seqTable, kmerTable, *k)
}

type Weight struct {
	weight float64
	index  int
}

type WeightList []Weight

func (self WeightList) Len() int {
	return len(self)
}

func (self *WeightList) Swap(i, j int) {
	(*self)[i], (*self)[j] = (*self)[j], (*self)[i]
}

func (self WeightList) Less(i, j int) bool {
	return self[i].weight > self[j].weight
}

func printFeature(out io.Writer, V, W, H matrix.Matrix, seqTable []string, kmerTable []kmerindex.Kmer, k int) {
	patternCount, seqCount := H.Dims()
	kmerCount, _ := W.Dims()

	hipats := make([]WeightList, seqCount)
	pats := make([]string, 0)

	for i := 0; i < patternCount; i++ {
		fmt.Fprintf(out, "Feature %v:\n", i)
		klist := make(WeightList, 0)
		for j := 0; j < kmerCount; j++ {
			klist = append(klist, Weight{weight: W.At(j, i), index: j})
		}
		sort.Sort(&klist)
		name := fmt.Sprint("[")
		for j := 0; j < len(klist); j++ {
			if klist[j].weight > 0 {
				ks, err := kmerindex.Format(kmerTable[klist[j].index], k, alphabet.DNA)
				if err != nil {
					panic(err)
				}
				name += fmt.Sprintf(" %s/%.3e ", ks, klist[j].weight)
			}
		}
		name += fmt.Sprint("]")
		pats = append(pats, name)
		fmt.Fprintln(out, name)

		slist := make(WeightList, 0)
		for j := 0; j < seqCount; j++ {
			slist = append(slist, Weight{weight: H.At(i, j), index: j})
			hipats[j] = append(hipats[j], Weight{weight: H.At(i, j), index: i})
		}

		sort.Sort(&slist)
		instances := ""
		for j := 0; j < len(slist); j++ {
			if slist[j].weight > 0 {
				instances += fmt.Sprintf("%s/%.3e\n", seqTable[slist[j].index], slist[j].weight)
			}
		}
		fmt.Fprintln(out, instances)

		fmt.Fprintln(out)
	}

	for j, seq := range hipats {
		sort.Sort(&seq)
		fmt.Fprintf(out, "%s/%err: %d\n", seqTable[j], seq[0].weight, seq[0].index)
	}
}
