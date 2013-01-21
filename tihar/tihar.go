package main

import (
	"code.google.com/p/biogo.matrix"
	"code.google.com/p/biogo/exp/seqio/fasta"
	"code.google.com/p/biogo/index/kmerindex"

	"code.google.com/p/biogo/exp/alphabet"
	"code.google.com/p/biogo/exp/seq/linear"
	"flag"
	"fmt"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"
	"sort"
	"time"
)

func main() {
	var (
		in                *fasta.Reader
		out, csv, profile *os.File
		err               error
	)

	inName := flag.String("in", "", "Filename for input to be factorised. Defaults to stdin.")
	outName := flag.String("out", "", "Filename for output. Defaults to stdout.")
	csvName := flag.String("csv", "", "Filename for csv output of feature details. Defaults to stdout.")
	k := flag.Int("k", 8, "kmer size to use.")
	cat := flag.Int("cat", 5, "number of categories.")
	iter := flag.Int("i", 1000, "iterations.")
	limit := flag.Int("time", 10, "time limit for NMF.")
	lo := flag.Int("lo", 1, "minimum number of kmer frequency to use in NMF.")
	hi := flag.Float64("hi", 0.9, "maximum proportion of kmer representation to use in NMF.")
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

	if *csvName == "" {
		fmt.Fprintln(os.Stderr, "Writing csv output to stdout.")
		csv = os.Stdout
	} else if csv, err = os.Create(*csvName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
	} else {
		fmt.Fprintf(os.Stderr, "Writing output to `%s'.\n", *csvName)
	}
	defer csv.Close()

	kmers := make(map[kmerindex.Kmer]int)
	positions := make(map[int]int)
	motifs := make(map[kmerindex.Kmer]map[int]map[string]bool)
	maxPos := 0

	for {
		if s, err := in.Read(); err != nil {
			break
		} else {
			if kindex, err := kmerindex.New(*k, s.(*linear.Seq)); err != nil {
				fmt.Fprintf(os.Stderr, "Error: %v.", err)
				os.Exit(1)
			} else {
				kindex.Build()
				index, _ := kindex.KmerIndex()
				for kmer, posList := range index {
					if _, ok := motifs[kmer]; !ok {
						motifs[kmer] = make(map[int]map[string]bool)
					}
					for _, pos := range posList {
						if _, ok := motifs[kmer][pos]; !ok {
							motifs[kmer][pos] = make(map[string]bool)
						}
						motifs[kmer][pos][string(s.Name())] = true
						kmers[kmer]++
						positions[pos]++
						if pos > maxPos {
							maxPos = pos
						}
					}
				}
			}
		}
	}

	kmerArray := make([][]float64, 0)
	kmerTable := make([]kmerindex.Kmer, 0)
	positionsTable := make(map[int]int)
	currPos := 0

	for kmer, count := range kmers {
		if count < *lo || float64(count)/float64(maxPos) > *hi {
			continue
		}
		row := make([]float64, currPos)
		for pos, seqs := range motifs[kmer] {
			if len(seqs) < *lo {
				continue
			}
			if i, ok := positionsTable[pos]; ok {
				row[i] += float64(len(motifs[kmer][pos]))
			} else {
				positionsTable[pos] = len(row)
				row = append(row, float64(len(motifs[kmer][pos])))
				currPos++
			}
		}
		kmerArray = append(kmerArray, row)
		kmerTable = append(kmerTable, kmerindex.Kmer(kmer))
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

	printFeature(out, csv, kMat, W, H, motifs, kmerTable, positionsTable, maxPos, *k)
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

func printFeature(out, csv *os.File, V, W, H matrix.Matrix, motifs map[kmerindex.Kmer]map[int]map[string]bool, kmerTable []kmerindex.Kmer, positionsTable map[int]int, maxPos, k int) {
	patternCount, posCount := H.Dims()
	kmerCount, _ := W.Dims()

	hipats := make([]WeightList, posCount)
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

		plist := make(WeightList, 0)
		for j := 0; j < posCount; j++ {
			plist = append(plist, Weight{weight: H.At(i, j), index: positionsTable[j]})
			hipats[j] = append(hipats[j], Weight{weight: H.At(i, j), index: i})
		}

		sort.Sort(&plist)
		instances := ""
		for j := 0; j < len(plist); j++ {
			if plist[j].weight > 0 {
				instances += fmt.Sprintf("%d/%.3e\n", plist[j].index, plist[j].weight)
			}
		}
		fmt.Fprintln(out, instances)

		fmt.Fprintln(out)
	}

	fmt.Fprint(csv, "position\tfeature\tweight")
	for j := 0; j < patternCount; j++ {
		fmt.Fprintf(csv, "\t%d", j)
	}
	fmt.Fprintln(csv)
	for i := 0; i <= maxPos; i++ {
		if pos, ok := positionsTable[i]; ok {
			all := ""
			for _, pat := range hipats[pos] {
				all += fmt.Sprintf("\t%err", pat.weight)
			}
			sort.Sort(&hipats[pos])
			if hipats[pos][0].weight > 0 {
				fmt.Fprintf(csv, "%d\t%d\t%err%s\n", i, hipats[pos][0].index, hipats[pos][0].weight, all)
			}
		}
	}
}
