// igor is a tool that takes pairwise alignment data as produced by PALS or krishna
// and generates repeat feature consensus sequences.
package main

import (
	"code.google.com/p/biogo.external/mafft"
	"code.google.com/p/biogo.external/muscle"
	"code.google.com/p/biogo.graph"
	"code.google.com/p/biogo/align/pals"
	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/featio/gff"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq"
	"code.google.com/p/biogo/seq/linear"
	"code.google.com/p/biogo/seq/multi"
	"code.google.com/p/biogo/seq/sequtils"

	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"net/http"
	_ "net/http/pprof"
	"os"
	"os/exec"
	"runtime"
	"runtime/pprof"
	"strings"
	"sync"
)

var bug debug

type debug bool

func (d debug) Printf(format string, v ...interface{}) (int, error) {
	if d {
		return fmt.Printf(format, v...)
	}
	return 0, nil
}

var aligner int

const (
	Muscle = iota + 1
	Mafft
)

type strandEdge struct {
	graph.Edge
	Strand seq.Strand
}

func main() {
	var (
		in  *gff.Reader
		out *gff.Writer
		err error
	)

	var (
		inName  = flag.String("in", "", "Filename for input.")
		outName = flag.String("out", "", "Filename for output. Defaults to stdout.")
		refName = flag.String("ref", "", "Filename of fasta file containing reference sequence.")

		epsilon = flag.Float64("epsilon", 0.15, "Tolerance for clustering.")
		effort  = flag.Int("effort", 5, "Number of attempts for clustering.")

		minFamily = flag.Int("famsize", 2, "Minimum number of clusters per family (must be >= 2).")

		cpuprofile = flag.String("cpuprofile", "", "Write cpu profile to this file.")
		webprofile = flag.String("webprofile", "", "Run web-based profiling on this host:port.")

		help = flag.Bool("help", false, "Print usage message.")
	)
	flag.IntVar(&aligner, "aligner", 1, "Which aligner to use: 1 - muscle, 2 - mafft")
	flag.BoolVar((*bool)(&bug), "debug", false, "Print graph generation information.")

	flag.Parse()
	if *minFamily < 2 {
		*minFamily = 2
	}

	if *help {
		flag.Usage()
		os.Exit(0)
	}

	if *webprofile != "" {
		go func() {
			log.Println(http.ListenAndServe(*webprofile, nil))
		}()
	}
	if *cpuprofile != "" {
		profile, err := os.Create(*cpuprofile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.", err)
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Writing CPU profile data to %s\n", *cpuprofile)
		pprof.StartCPUProfile(profile)
		defer pprof.StopCPUProfile()
	}

	if *inName == "" {
		fmt.Fprintln(os.Stderr, "Reading PALS features from stdin.")
		in = gff.NewReader(os.Stdin)
	} else if f, err := os.Open(*inName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		defer f.Close()
		in = gff.NewReader(f)
		fmt.Fprintf(os.Stderr, "Reading PALS features from `%s'.\n", *inName)
	}

	if *outName == "" {
		fmt.Fprintln(os.Stderr, "Writing to stdout.")
		out = gff.NewWriter(os.Stdout, 60, false)
		out.Precision = 2
	} else if f, err := os.Create(*outName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		defer f.Close()
		out = gff.NewWriter(f, 60, true)
		out.Precision = 2
		fmt.Fprintf(os.Stderr, "Writing to `%s'.\n", *outName)
	}
	defer out.Close()

	fmt.Fprintf(os.Stderr, "Building data structures.\n")
	fmt.Fprintf(os.Stderr, " Generating piles ...\n")
	piler := pals.NewPiler(0)
	for {
		rep, err := in.Read()
		if err != nil {
			if err != io.EOF {
				fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
				os.Exit(1)
			}
			break
		}

		p, err := pals.ExpandFeature(rep.(*gff.Feature))
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
			os.Exit(1)
		}
		piler.Add(p)

	}
	piles, err := piler.Piles(nil)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
		os.Exit(1)
	}

	fmt.Fprintf(os.Stderr, " Subclustering piles ...")
	cls := make([]*Cluster, len(piles))
	wg := &sync.WaitGroup{}
	q := make(chan struct{}, runtime.GOMAXPROCS(0))
	for i, p := range piles {
		wg.Add(1)
		q <- struct{}{}
		go func(i int, p *pals.Pile) {
			defer func() { <-q; wg.Done() }()
			cls[i] = BuildCluster(p, *epsilon, *effort)
		}(i, p)
	}
	wg.Wait()
	var p int
	for _, cl := range cls {
		p += len(cl.Piles)
	}
	fmt.Fprintf(os.Stderr, " generated %d subpiles from %d piles.\n", p, len(cls))

	fmt.Fprintf(os.Stderr, " Connected components ...\n")
	g := graph.NewUndirected()

	type feature struct {
		name     string
		from, to int
	}
	added := map[feature]graph.Node{}

	for _, cl := range cls {
		bug.Printf("\n%v %v\n", cl.Name(), cl.Description())
		for _, pi := range cl.Piles {
			if pi.Node == nil {
				n, ok := added[feature{pi.Name(), pi.Start(), pi.End()}]
				if !ok {
					n = g.NewNode()
					pi.Node = n
					g.Add(pi)
					added[feature{pi.Name(), pi.Start(), pi.End()}] = n
				} else {
					pi.Node = n
				}
			}
			bug.Printf("\t%v %v %v\n",
				pi.Name(), pi.Description(), pi.Node)
			for _, im := range pi.Images {
				impp := im.B.Location().(*pals.Pile)
				if impp.Node == nil {
					n, ok := added[feature{impp.Name(), impp.Start(), impp.End()}]
					if !ok {
						n = g.NewNode()
						impp.Node = n
						g.Add(impp)
						added[feature{pi.Name(), pi.Start(), pi.End()}] = n
					} else {
						impp.Node = n
					}
				}
				if pi.Node != impp.Node {
					// TODO Only make non-redundant connections:
					//  near self alignments should not be connected.
					e := strandEdge{Edge: graph.NewEdge(), Strand: im.Strand}
					e.SetWeight(float64(im.Score))
					err := g.ConnectWith(pi, impp, e)
					if err != nil {
						panic(fmt.Sprintf("internal error: %v", err))
					}
				}
				bug.Printf("\t\t%v %v %v %v\n",
					impp.Name(), impp.Description(), im.Score, impp.Node)
			}
		}
	}

	cc := g.ConnectedComponents(graph.EdgeFilter(func(e graph.Edge) bool {
		h, t := e.Head().(*pals.Pile), e.Tail().(*pals.Pile)
		switch {
		case h.Strand == 0 && t.Strand == 0:
			h.Strand = 1
			t.Strand = e.(strandEdge).Strand
		case t.Strand == 0:
			t.Strand = h.Strand * e.(strandEdge).Strand
		case h.Strand == 0:
			h.Strand = t.Strand * e.(strandEdge).Strand
		default:
			if s := e.(strandEdge).Strand; h.Strand != t.Strand*s {
				if bug {
					fmt.Printf("igor: inconsistency in strands: %d != %d under T(%d)\n", h.Strand, t.Strand, s)
				} else {
					panic(fmt.Sprintf("igor: inconsistency in strands: %d != %d under T(%d)\n", h.Strand, t.Strand, s))
				}
			}
		}
		return true
	}))

	if *refName == "" {
		featureList(cc, *minFamily)
	} else {
		sequences(cc, *minFamily, *refName)
	}
}

func featureList(cc []graph.Nodes, minFamily int) {
	for fi, fam := range cc {
		if len(fam) < minFamily {
			continue
		}
		fmt.Printf("Family %d (%d members):\n", fi, len(fam))
		for _, sfam := range fam {
			var (
				sfp    = sfam.(*pals.Pile)
				cl     = sfp.Location()
				contig = cl.Location()
			)
			fmt.Printf(" %v %s %d %d %d\n", contig.Name(), sfp.Description(), sfp.Start(), sfp.End(), sfp.Strand)
			for _, pi := range sfp.Images {
				var (
					pl     = pi.B.Location()
					cl     = pl.Location()
					contig = cl.Location()
				)
				fmt.Printf(" - %s %s %d %d %v\n", contig.Name(), pi.B.Description(), pi.B.Start(), pi.B.End(), pi.Strand)
			}
		}
		fmt.Println()
	}
}

func sequences(cc []graph.Nodes, minFamily int, refName string) {
	refStore := map[string]*linear.Seq{}
	f, err := os.Open(refName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	}
	r := fasta.NewReader(f, &linear.Seq{Annotation: seq.Annotation{Alpha: alphabet.DNA}})
	for {
		s, err := r.Read()
		if err != nil {
			break
		}
		refStore[s.Name()] = s.(*linear.Seq)
	}

	var rc = make(chan struct {
		fam  string
		cons *linear.QSeq
	})
	go func() {
		var (
			wg = &sync.WaitGroup{}
			q  = make(chan struct{}, runtime.GOMAXPROCS(0))
		)
		for fi, fam := range cc {
			if len(fam) < minFamily {
				continue
			}
			wg.Add(1)
			q <- struct{}{}
			go func(fi int, fam graph.Nodes) {
				defer func() { <-q; wg.Done() }()
				var s string
				for _, sfam := range fam {
					var (
						sfp    = sfam.(*pals.Pile)
						cl     = sfp.Location()
						contig = cl.Location()
					)
					ss := *refStore[contig.Name()]
					sequtils.Truncate(&ss, refStore[contig.Name()], sfp.Start(), sfp.End())
					if sfp.Strand < 0 {
						ss.RevComp()
					}
					ss.ID = contig.Name()
					ss.Desc = fmt.Sprintf("%s in %s (%d %d %v)",
						sfp.Description(), cl.Name(), sfp.Start(), sfp.End(), sfp.Strand)
					s = fmt.Sprintf("%s%60a\n", s, &ss)
				}
				var m *exec.Cmd
				switch aligner {
				case Muscle:
					m, err = muscle.Muscle{Quiet: true}.BuildCommand()
				case Mafft:
					m, err = mafft.Mafft{InFile: "-", Auto: true, Quiet: true}.BuildCommand()
				default:
					panic("igor: no valid aligner specified")
				}
				if err != nil {
					panic(err)
				}
				m.Stdin = strings.NewReader(s)
				m.Stdout = &bytes.Buffer{}
				m.Run()
				c := genConsensus(m.Stdout.(io.Reader), fmt.Sprintf("Family_%d", fi))
				c.Desc = fmt.Sprintf("(%d members)", len(fam))
				rc <- struct {
					fam  string
					cons *linear.QSeq
				}{
					fam:  fmt.Sprintf("Family_%d (%d members):\n%s", fi, len(fam), m.Stdout),
					cons: c,
				}
			}(fi, fam)
		}
		wg.Wait()
		close(rc)
	}()
	for r := range rc {
		// fmt.Print(r.fam)
		r.cons.Threshold = seq.DefaultQphred
		r.cons.QFilter = seq.CaseFilter
		fmt.Printf("%60a\n", r.cons)
	}
}

func genConsensus(in io.Reader, id string) *linear.QSeq {
	var (
		r = fasta.NewReader(in, &linear.Seq{
			Annotation: seq.Annotation{Alpha: alphabet.DNA},
		})

		ms = &multi.Multi{
			Annotation:     seq.Annotation{ID: id},
			ColumnConsense: seq.DefaultQConsensus,
		}
	)
	for {
		s, err := r.Read()
		if err != nil {
			break
		}
		ms.Add(s)
	}
	return ms.Consensus(true)
}
