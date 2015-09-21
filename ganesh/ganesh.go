package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/gonum/matrix/mat64"

	"github.com/kortschak/nmf"
)

func main() {
	var (
		in      *bufio.Reader
		out     *bufio.Writer
		profile *os.File
		e       error
	)

	inName := flag.String("in", "", "Filename for input to be factorised. Defaults to stdin.")
	outName := flag.String("out", "", "Filename for output. Defaults to stdout.")
	transpose := flag.Bool("t", false, "Transpose columns and rows.")
	sep := flag.String("sep", "\t", "Column delimiter.")
	cat := flag.Int("cat", 5, "number of categories.")
	iter := flag.Int("i", 1000, "iterations.")
	rep := flag.Int("rep", 1, "Resample replicates.")
	limit := flag.Duration("time", 10*time.Second, "time limit for NMF.")
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
		if profile, e = os.Create(*cpuprofile); e != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.", e)
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Writing CPU profile data to %s\n", *cpuprofile)
		pprof.StartCPUProfile(profile)
		defer pprof.StopCPUProfile()
	}

	if *inName == "" {
		fmt.Fprintln(os.Stderr, "Reading table from stdin.")
		in = bufio.NewReader(os.Stdin)
	} else if f, err := os.Open(*inName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		defer f.Close()
		in = bufio.NewReader(f)
		fmt.Fprintf(os.Stderr, "Reading table from `%s'.\n", *inName)
	}

	if *outName == "" {
		fmt.Fprintln(os.Stderr, "Writing output to stdout.")
		out = bufio.NewWriter(os.Stdout)
	} else if f, err := os.Create(*outName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		defer f.Close()
		out = bufio.NewWriter(f)
		fmt.Fprintf(os.Stderr, "Writing output to `%s'.\n", *outName)
	}
	defer out.Flush()

	var (
		colNames, rowNames []string

		r    int
		data []float64
	)
	if line, err := in.ReadString('\n'); err != nil {
		fmt.Fprintln(os.Stderr, "No table to read\n")
		os.Exit(1)
	} else {
		line = strings.TrimSpace(line)
		colNames = strings.Split(line, "\t")
		colNames = colNames[1:]
	}

	for count := 1; ; count++ {
		if line, err := in.ReadString('\n'); err != nil {
			break
		} else {
			line = strings.TrimSpace(line)
			if row := strings.Split(line, *sep); len(row) != len(colNames)+1 {
				fmt.Fprintf(os.Stderr, "Table row mismatch at line %d.\n", count)
				os.Exit(1)
			} else {
				rowData := make([]float64, len(row)-1)
				for i, val := range row[1:] {
					if rowData[i], e = strconv.ParseFloat(val, 64); e != nil {
						fmt.Fprintf(os.Stderr, "Float conversion error %v at line %d element %d.\n", e, count, i)
						os.Exit(1)
					}
				}
				rowNames = append(rowNames, row[0])
				data = append(data, rowData...)
				r++
			}
		}
	}

	mat := mat64.NewDense(r, len(colNames), data)

	var nonZero float64
	f := func(_, _ int, v float64) float64 {
		if v != 0 {
			nonZero++
		}
		return v
	}
	mat.Apply(f, mat)

	if *transpose {
		colNames, rowNames = rowNames, colNames
		mat.Clone(mat.T())
	}
	r, c := mat.Dims()

	density := nonZero / float64(r*c)

	if *seed == -1 {
		*seed = time.Now().UnixNano()
	}
	fmt.Fprintf(os.Stderr, "Using %v as random seed.\n", *seed)
	rand.Seed(*seed)

	fmt.Fprintf(os.Stderr, "Dimensions of matrix = (%v, %v)\nDensity = %.3f %%\n%v\n", r, c, (density)*100, mat)

	for run := 0; run < *rep; run++ {
		if *rep > 1 {
			fmt.Fprintf(os.Stderr, "Replicate #%d\n", run+1)
		}

		posNorm := func(_, _ int, _ float64) float64 { return math.Abs(rand.NormFloat64()) }

		Wo := mat64.NewDense(r, *cat, nil)
		Wo.Apply(posNorm, Wo)

		Ho := mat64.NewDense(*cat, c, nil)
		Ho.Apply(posNorm, Ho)

		W, H, ok := nmf.Factors(mat, Wo, Ho, nmf.Config{Tolerance: *tol, MaxIter: *iter, Limit: *limit})

		fmt.Fprintf(os.Stderr, "norm(H) = %v norm(W) = %v\n\nFinished = %v\n\n", H.Norm(0), W.Norm(0), ok)

		printFeature(out, run, mat, W, H, rowNames, colNames)
	}
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

func printFeature(out io.Writer, run int, V, W, H *mat64.Dense, rowNames, colNames []string) {
	patternCount, colCount := H.Dims()
	rowCount, _ := W.Dims()

	hipats := make([]WeightList, colCount)
	pats := make([]string, 0)

	for i := 0; i < patternCount; i++ {
		rlist := make(WeightList, 0)
		for j := 0; j < rowCount; j++ {
			rlist = append(rlist, Weight{weight: W.At(j, i), index: j})
		}
		sort.Sort(&rlist)
		name := []string{}
		for j := 0; j < len(rlist); j++ {
			if rlist[j].weight > 0 {
				name = append(name, fmt.Sprintf("%s/%.3e", rowNames[rlist[j].index], rlist[j].weight))
			}
		}
		nameString := strings.Join(name, ",")
		pats = append(pats, nameString)

		clist := make(WeightList, 0)
		for j := 0; j < colCount; j++ {
			clist = append(clist, Weight{weight: H.At(i, j), index: j})
			hipats[j] = append(hipats[j], Weight{weight: H.At(i, j), index: i})
		}

		sort.Sort(&clist)
		instances := []string{}
		for j := 0; j < len(clist); j++ {
			if clist[j].weight > 0 {
				instances = append(instances, fmt.Sprintf("%s/%.3e", colNames[clist[j].index], clist[j].weight))
			}
		}
		instanceString := strings.Join(instances, ",")
		fmt.Fprintf(out, "%d\t[%s]\t(%s)\n", run, nameString, instanceString)
	}

	for j, col := range hipats {
		sort.Sort(&col)
		fmt.Fprintf(os.Stderr, "%s/%e: %d\n", colNames[j], col[0].weight, col[0].index)
	}
}
