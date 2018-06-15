// krishna is a pure Go implementation of Edgar and Myers PALS tool.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"net/http"
	_ "net/http/pprof"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/biogo/biogo/align/pals"
	"github.com/biogo/biogo/align/pals/filter"
	"github.com/biogo/biogo/morass"
)

const (
	timeFormat = "20060102150405-Mon"
)

var (
	pid           = os.Getpid()
	mem           *uintptr
	profile       *os.File
	queryName     string
	targetName    string
	selfCompare   bool
	sameStrand    bool
	outFile       string
	trapFile      bool
	maxK          int
	minHitLen     int
	minId         float64
	dpMinHitLen   int
	dpMinId       float64
	tubeOffset    int
	tmpDir        string
	tmpChunk      int
	tmpConcurrent bool
	threads       int
	maxMem        uint64
	logToFile     bool
	debug         bool
	verbose       bool
	cpuprofile    string
	webprofile    string
	logger        *log.Logger
)

func init() {
	flag.StringVar(&queryName, "query", "", "Filename for query sequence.")
	flag.StringVar(&targetName, "target", "", "Filename for target sequence.")
	flag.BoolVar(&selfCompare, "self", false, "Is this a self comparison?")
	flag.BoolVar(&sameStrand, "same", false, "Only compare same strand")

	flag.StringVar(&outFile, "out", "", "File to send output to.")
	flag.BoolVar(&trapFile, "traps", false, "Specifies whether to keep trapezoid seeds.")

	flag.IntVar(&maxK, "k", -1, "Maximum kmer length (negative indicates automatic detection based on architecture).")
	flag.IntVar(&minHitLen, "filtlen", 400, "Minimum hit length for filter.")
	flag.Float64Var(&minId, "filtid", 0.94, "Minimum hit identity for filter.")
	flag.IntVar(&dpMinHitLen, "dplen", 0, "Minimum hit length for aligner.")
	flag.Float64Var(&dpMinId, "dpid", 0, "Minimum hit identity for aligner.")
	flag.IntVar(&tubeOffset, "tubeoffset", 0, "Tube offset - 0 indicate autotune.")

	flag.StringVar(&tmpDir, "tmp", "", "Path for temporary files.")
	flag.IntVar(&tmpChunk, "chunk", 100e6, "Chunk size for morass.")
	flag.BoolVar(&tmpConcurrent, "tmpcon", false, "Process morass concurrently.")

	flag.IntVar(&threads, "threads", 1, "Number of threads to use for alignment.")
	flag.Uint64Var(&maxMem, "mem", 0, "Maximum nominal memory - 0 indicates unlimited.")

	flag.BoolVar(&logToFile, "log", false, "Log to file.")
	flag.BoolVar(&debug, "debug", false, "Include file names/lines in log.")
	flag.BoolVar(&verbose, "v", false, "Log additional information.")

	flag.StringVar(&cpuprofile, "cpuprofile", "", "write cpu profile to this file.")
	flag.StringVar(&webprofile, "webprofile", "", "Run web-based profiling on this host:port.")

	help := flag.Bool("help", false, "Print this help message.")

	flag.Parse()

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	if maxMem == 0 {
		mem = nil
	} else {
		*mem = uintptr(maxMem)
	}

	if threads > runtime.GOMAXPROCS(0) {
		runtime.GOMAXPROCS(threads)
	}
}

func initLog(fileName string) {
	var w io.Writer = os.Stderr
	if fileName != "" {
		file, err := os.Create(fileName)
		if err == nil {
			fmt.Fprintln(file, strings.Join(os.Args, " "))
			w = io.MultiWriter(os.Stderr, file)
		} else {
			fmt.Printf("Error: Could not open log file: %v", err)
			os.Exit(1)
		}
	}

	logger = log.New(w, fmt.Sprintf("%s:", filepath.Base(os.Args[0])), log.Flags())
	if debug {
		logger.SetFlags(log.Flags() | log.Lshortfile)
	}
}

func main() {
	if webprofile != "" {
		go func() {
			log.Println(http.ListenAndServe(webprofile, nil))
		}()
	}
	if cpuprofile != "" {
		profile, err := os.Create(cpuprofile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.", err)
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Writing CPU profile data to %s\n", cpuprofile)
		pprof.StartCPUProfile(profile)
		defer pprof.StopCPUProfile()
	}

	if logToFile {
		initLog("krishna-" + time.Now().Format(timeFormat) + "-" + strconv.Itoa(pid) + ".log")
	} else {
		initLog("")
	}

	logger.Println(os.Args)
	var target, query *pals.Packed
	if targetName != "" {
		var err error
		target, err = packSequence(targetName)
		if err != nil {
			log.Fatalf("Internal error: %v", err)
		}
		if target.Len() == 0 {
			log.Fatal("Target sequence is zero length.")
		}
	} else {
		logger.Fatalln("No target provided.")
	}

	if !selfCompare {
		if queryName != "" {
			var err error
			query, err = packSequence(queryName)
			if err != nil {
				log.Fatalf("Internal error: %v", err)
			}
			if query.Len() == 0 {
				log.Fatal("Query sequence is zero length.")
			}
		} else {
			logger.Fatalln("No query provided in non-self comparison.")
		}
	} else {
		query = target
	}

	var writer *pals.Writer
	if outFile == "" {
		writer = pals.NewWriter(os.Stdout, 2, 60, false)
	} else {
		out, err := os.Create(outFile)
		if err != nil {
			log.Fatalf("Could not open output file: %v", err)
		}
		defer out.Close()
		buf := bufio.NewWriter(out)
		defer buf.Flush()
		writer = pals.NewWriter(buf, 2, 60, false)
	}

	if maxK > 0 {
		pals.MaxKmerLen = maxK
	}

	m := func() *morass.Morass {
		m, err := morass.New(filter.Hit{}, "krishna_"+strconv.Itoa(pid)+"_", tmpDir, tmpChunk, tmpConcurrent)
		if err != nil {
			logger.Fatalf("Error: %v", err)
		}
		return m
	}
	pa := []*pals.PALS{pals.New(target.Seq, query.Seq, selfCompare, m(), tubeOffset, mem, logger)}
	if threads > 1 {
		pa = append(pa, pals.New(target.Seq, query.Seq, selfCompare, m(), tubeOffset, mem, logger))
	}

	if err := pa[0].Optimise(minHitLen, minId); err != nil {
		logger.Fatalf("Error: %v", err)
	}
	if dpMinHitLen != 0 {
		pa[0].DPParams.MinHitLength = dpMinHitLen
	}
	if dpMinId != 0 {
		pa[0].DPParams.MinId = dpMinId
	}

	logger.Printf("Using filter parameters:")
	logger.Printf("\tWordSize = %d", pa[0].FilterParams.WordSize)
	logger.Printf("\tMinMatch = %d", pa[0].FilterParams.MinMatch)
	logger.Printf("\tMaxError = %d", pa[0].FilterParams.MaxError)
	logger.Printf("\tTubeOffset = %d", pa[0].FilterParams.TubeOffset)
	logger.Printf("\tAvg List Length = %.3f", pa[0].AvgIndexListLength(pa[0].FilterParams))
	logger.Printf("Using dynamic programming parameters:")
	logger.Printf("\tMinLen = %d", pa[0].DPParams.MinHitLength)
	logger.Printf("\tMinID = %.1f%%", pa[0].DPParams.MinId*100)
	logger.Printf("Estimated minimum memory required = %dMiB", pa[0].MemRequired(pa[0].FilterParams)/(1<<20))
	logger.Printf("Building index for %s", target.ID)

	if err := pa[0].BuildIndex(); err != nil {
		logger.Fatalf("Error: %v", err)
	}
	if threads > 1 {
		pa[1].Share(pa[0])
	}

	both := !sameStrand
	wg := &sync.WaitGroup{}
	for i, comp := range [...]bool{false, true} {
		if threads > 1 && both {
			wg.Add(1)
			go func(p *pals.PALS, comp bool) {
				defer wg.Done()
				hits, err := p.Align(comp)
				if err != nil {
					logger.Fatalf("Error: %v", err)
				}
				if trapFile {
					logger.Println("Writing trapezoid data")
					err = WriteTraps(comp, p.Trapezoids())
					if err != nil {
						logger.Fatalf("Error: %v", err)
					}
				}

				logger.Println("Writing results")
				n, err := WriteDPHits(writer, target, query, hits, comp)
				if err != nil {
					logger.Fatalf("Error: %v.", err)
				}
				logger.Printf("Wrote hits (%v bytes)", n)
			}(pa[i], comp)
		} else {
			if comp {
				logger.Println("Working on complementary strands")
			} else {
				logger.Println("Working on self strand")
			}
			if both || !comp {
				hits, err := pa[0].Align(comp)
				if err != nil {
					logger.Fatalf("Error: %v", err)
				}
				if trapFile {
					logger.Println("Writing trapezoid data")
					err = WriteTraps(comp, pa[0].Trapezoids())
					if err != nil {
						logger.Fatalf("Error: %v", err)
					}
				}

				logger.Println("Writing results")
				n, err := WriteDPHits(writer, target, query, hits, comp)
				if err != nil {
					logger.Fatalf("Error: %v.", err)
				}
				logger.Printf("Wrote hits (%v bytes)", n)
			}
		}
	}
	wg.Wait()

	for _, p := range pa {
		p.CleanUp()
	}

	logger.Print("Finished.")
}
