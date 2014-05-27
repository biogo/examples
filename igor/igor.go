// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// igor is a tool that takes pairwise alignment data as produced by PALS or krishna
// and reports repeat feature family groupings in JSON format.
package main

import (
	"code.google.com/p/biogo.examples/igor/igor"

	"code.google.com/p/biogo.graph"
	"code.google.com/p/biogo/align/pals"
	"code.google.com/p/biogo/io/featio/gff"
	"code.google.com/p/biogo/seq"

	"bufio"
	"encoding/json"
	"flag"
	"io"
	"log"
	"os"
)

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func within(alpha float64, short, long int) bool {
	return float64(short) >= float64(long)*(1-alpha)
}

var (
	inName  string
	outName string

	landscapeDir string

	classic bool

	pileDiff  float64
	imageDiff float64

	band          float64
	mergeOverlap  int
	removeOverlap float64
	requiredCover float64
	strictness    int

	threads int
	bug     bool
)

func init() {
	flag.StringVar(&inName, "in", "", "Filename for input.")
	flag.StringVar(&outName, "out", "", "Filename for output. Defaults to stdout.")
	flag.StringVar(&landscapeDir, "landscapes", "", "Directory to output landscape data (deletes existing directory).")

	flag.Float64Var(&band, "band", 0.05, "Kernel bandwidth as fraction of pile length.")
	flag.Float64Var(&pileDiff, "pile-diff", 0.05, "Fractional length difference tolerance between piles.")
	flag.Float64Var(&imageDiff, "image-diff", 0.05, "Fractional length difference tolerance for images and piles.")
	flag.IntVar(&mergeOverlap, "merge-overlap", 0, "Threshold for merging adjacent images into piles.")
	flag.Float64Var(&removeOverlap, "remove-overlap", 0.95, "Fractional pile overlap threshold for removal in clustering.")
	flag.Float64Var(&requiredCover, "target-coverage", 0.95, "Fractional pile coverage threshold.")
	flag.IntVar(&strictness, "overlap-strictness", 0, "Keep overlapping sub piles (0-2).")

	flag.IntVar(&threads, "threads", 4, "Number of parallel clustering threads to use.")
	flag.BoolVar(&bug, "debug", false, "Print graph generation information.")

	flag.BoolVar(&classic, "classic", false, "Run a reasonable approximation of the C implementation of PILER.")

	help := flag.Bool("help", false, "Print usage message.")

	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	if strictness < 0 || strictness > 2 {
		flag.Usage()
		os.Exit(1)
	}

	if landscapeDir != "" {
		fi, err := os.Stat(landscapeDir)
		if err != nil && !os.IsNotExist(err) {
			log.Fatalf("failed to stat landscape directory: %v", err)
		}
		if fi != nil && !fi.IsDir() {
			log.Fatalf("will not delete non-directory file %q", landscapeDir)
		}
		err = os.RemoveAll(landscapeDir)
		if err != nil {
			log.Fatalf("failed to delete landscape directory %q: %v", landscapeDir, err)
		}
		err = os.Mkdir(landscapeDir, 0755)
		if err != nil {
			log.Fatalf("failed to create landscape directory %q: %v", landscapeDir, err)
		}
	}
}

func main() {
	var (
		in  *gff.Reader
		out io.Writer
	)

	if inName == "" {
		log.Println("reading PALS features from stdin")
		in = gff.NewReader(os.Stdin)
	} else if f, err := os.Open(inName); err != nil {
		log.Fatalf("error: %v", err)
	} else {
		defer f.Close()
		in = gff.NewReader(f)
		log.Printf("reading PALS features from %q\n", inName)
	}

	if outName == "" {
		log.Println("writing to stdout")
		out = os.Stdout
	} else if f, err := os.Create(outName); err != nil {
		log.Fatalf("error: %v", err)
	} else {
		defer f.Close()
		buf := bufio.NewWriter(f)
		defer buf.Flush()
		out = buf
		log.Printf("writing to %q\n", outName)
	}

	log.Println("generating piles ... piling.")
	var pf pals.PairFilter
	if classic {
		pf = func(p *pals.Pair) bool {
			// Test whether each image's length is within pileDiff of its pile's length.
			imageOK := within(imageDiff, p.A.Len(), p.A.Loc.Len()) && within(imageDiff, p.B.Len(), p.B.Loc.Len())

			// Test whether the shorter pile's is within pileDiff of the longer pile's length.
			pilesOK := within(pileDiff, min(p.A.Loc.Len(), p.B.Loc.Len()), max(p.A.Loc.Len(), p.B.Loc.Len()))

			return imageOK && pilesOK
		}
	}
	piles, err := igor.Piles(in, mergeOverlap, pf)
	if err != nil {
		log.Fatalf("piling error:", err)
	}

	var clusters [][]*pals.Pile
	if classic {
		clusters = [][]*pals.Pile{piles}
		log.Printf("generated %d piles.\n", len(piles))
	} else {
		log.Printf("subclustering piles (%d to do) ...\n", len(piles))
		var n int
		n, clusters = igor.Cluster(piles, igor.ClusterConfig{
			BandWidth:         band,
			RequiredCover:     requiredCover,
			OverlapStrictness: byte(strictness),
			OverlapThresh:     removeOverlap,
			LandscapeDir:      landscapeDir,
			Threads:           threads,
		})
		log.Printf("generated %d subpiles.\n", n)
	}

	log.Println("finding connected components ...")
	cc := igor.Group(clusters, igor.GroupConfig{
		pileDiff,
		imageDiff,
		classic,
		bug,
	})
	log.Printf("%d remaining connected components\n", len(cc))

	err = writeJSON(cc, out)
	if err != nil {
		log.Fatalf("error: %v", err)
	}
}

func writeJSON(cc []graph.Nodes, w io.Writer) error {
	type feat struct {
		C string
		S int
		E int
		O seq.Strand
	}
	var (
		a feat
		f []feat
		j = json.NewEncoder(w)
	)

	seen := make(map[feat]struct{})
	var fi int
	for _, fam := range cc {
		for _, p := range fam {
			pile := p.(*pals.Pile)

			a.C = pile.Location().Name()
			a.S = pile.Start()
			a.E = pile.End()
			a.O = pile.Strand
			if _, ok := seen[a]; !ok && pile.Loc != nil {
				seen[a] = struct{}{}
				f = append(f, a)
			}

			for _, im := range pile.Images {
				partner := im.Mate().Location().(*pals.Pile)
				if partner.Loc == nil || partner.Strand == seq.None && len(f) > 0 {
					continue
				}

				a.C = partner.Location().Name()
				a.S = partner.Start()
				a.E = partner.End()
				a.O = partner.Strand
				if _, ok := seen[a]; ok {
					continue
				}

				seen[a] = struct{}{}
				f = append(f, a)
			}
		}
		switch len(f) {
		case 0:
			continue
		case 1:
			f[0].O = seq.Plus
		default:
			for i := 0; i < len(f); {
				if f[i].O == seq.None {
					f[i], f = f[len(f)-1], f[:len(f)-1]
				} else {
					i++
				}
			}
		}
		if len(f) < 2 {
			continue
		}
		log.Printf("Family#%d (%d members)\n", fi, len(f))
		err := j.Encode(f)
		if err != nil {
			return err
		}
		f = f[:0]
		fi++
	}

	return nil
}
