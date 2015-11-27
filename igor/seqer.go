// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

// seqer returns multiple fasta sequences corresponding to feature intervals
// described in the JSON output from igor. It will also produce fastq consensus
// sequence output from one of MUSCLE or MAFFT.
package main

import (
	"bytes"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/seq/multi"
	"github.com/biogo/biogo/seq/sequtils"
	"github.com/biogo/external/mafft"
	"github.com/biogo/external/muscle"
)

type feat struct {
	Chr    string     `json:"C"`
	Start  int        `json:"S"`
	End    int        `json:"E"`
	Orient seq.Strand `json:"O"`
}

var (
	refName    string
	dir        string
	aligner    string
	maxFam     int
	subSample  bool
	minFamily  int
	lengthFrac float64
	threads    int
	consFasta  bool
)

func main() {
	flag.IntVar(&maxFam, "maxFam", 0, "maxFam indicates maximum family size considered (0 == no limit).")
	flag.BoolVar(&subSample, "subsample", false, "Choose maxFam members of a family if the family has more than maxFam members.")
	flag.IntVar(&minFamily, "famsize", 2, "Minimum number of clusters per family (must be >= 2).")
	flag.IntVar(&threads, "threads", 1, "Number of concurrent aligner instances to run.")
	flag.StringVar(&refName, "ref", "", "Filename of fasta file containing reference sequence.")
	flag.StringVar(&aligner, "aligner", "", "Aligner to use to generate consensus (muscle or mafft).")
	flag.BoolVar(&consFasta, "fasta", false, "Output consensus as fasta with quality case filtering.")
	flag.Float64Var(&lengthFrac, "minLen", 0, "Minimum proportion of longest family member.")
	flag.StringVar(&dir, "dir", "", "Target directory for output. If not empty dir is deleted first.")
	flag.Parse()

	if len(flag.Args()) != 1 {
		fmt.Fprintln(os.Stderr, "Need single input gff file (output from gffer or victor).")
		flag.Usage()
		os.Exit(1)
	}
	if refName == "" {
		fmt.Fprintln(os.Stderr, "Need reference.")
		flag.Usage()
	}
	if minFamily < 2 {
		minFamily = 2
	}

	if threads < 1 {
		threads = 1
	}
	manager.limit = make(chan struct{}, threads)

	if dir != "" {
		err := os.RemoveAll(dir)
		if err != nil {
			log.Fatalf("failed to remove target directory: %v", err)
		}
		err = os.Mkdir(dir, os.ModeDir|0750)
		if err != nil {
			log.Fatalf("failed to create target directory: %v", err)
		}
	}

	refStore := getReference(refName)

	f, err := os.Open(flag.Args()[0])
	if err != nil {
		log.Printf("error: could not open %s to read %v", flag.Args()[0], err)
	}
	defer f.Close()

	var v []*gff.Feature
	r := familyReader{r: gff.NewReader(f)}
	for {
		err := r.readInto(&v)
		if err != nil {
			if err != io.EOF {
				log.Fatalf("failed reading input file %q: %v", flag.Args()[0], err)
			}
			break
		}

		if len(v) < minFamily {
			continue
		}

		fam, err := strconv.Atoi(v[0].FeatAttributes.Get("Family"))
		if err != nil {
			log.Fatalf("failed to parse family id %q: %v", v[0].FeatAttributes, err)
		}

		var maxLen int
		for _, f := range v {
			if l := f.Len(); l > maxLen {
				maxLen = l
			}
		}
		lenThresh := int(float64(maxLen) * lengthFrac)

		var validLengthed int
		for _, f := range v {
			if f.Len() >= lenThresh {
				validLengthed++
			}
		}
		if maxFam != 0 && !subSample && validLengthed > maxFam {
			continue
		}

		var out *os.File
		if dir != "" {
			file := fmt.Sprintf("family%06d.mfa", fam)
			out, err = os.Create(filepath.Join(dir, file))
			if err != nil {
				log.Fatalf("failed to create %s: %v", file, err)
			}
		}

		if subSample {
			// Shuffle first.
			w := make([]*gff.Feature, 0, len(v))
			for _, j := range rand.Perm(len(v)) {
				w = append(w, v[j])
			}
			v = w
		}

		var sampled int
		for id, f := range v {
			if f.Len() < lenThresh {
				continue
			}
			if sampled++; subSample && sampled > maxFam {
				break
			}
			ss := *refStore[f.SeqName]
			sequtils.Truncate(&ss, refStore[f.SeqName], f.FeatStart, f.FeatEnd)
			if f.FeatStrand == seq.Minus {
				ss.RevComp()
			}
			ss.ID = fmt.Sprintf("family%06d_member%04d", fam, id)
			ss.Desc = fmt.Sprintf("%s:%d-%d %v (%d members - %d members within %.2f of maximum length) %v",
				f.SeqName, f.FeatStart, f.FeatEnd, f.FeatStrand, len(v), validLengthed, lengthFrac, f.FeatAttributes,
			)
			if dir == "" {
				fmt.Printf("%60a\n", &ss)
			} else {
				fmt.Fprintf(out, "%60a\n", &ss)
			}
		}
		if dir == "" {
			fmt.Println()
		} else {
			file := out.Name()
			out.Close()
			fam, lv, validLengthed, lengthFrac := fam, len(v), validLengthed, lengthFrac
			acquire()
			go func() {
				defer release()
				if aligner != "" {
					c, err := consensus(file, aligner)
					if err != nil {
						log.Printf("failed to generate consensus for family%06d: %v", fam, err)
					} else {
						c.ID = fmt.Sprintf("family%06d_consensus", fam)
						c.Desc = fmt.Sprintf("(%d members - %d members within %.2f of maximum length)",
							lv, validLengthed, lengthFrac,
						)
						c.Threshold = 42
						c.QFilter = seq.CaseFilter
						file := fmt.Sprintf("family%06d_consensus.fq", fam)
						out, err := os.Create(filepath.Join(dir, file))
						if err != nil {
							log.Printf("failed to create %s: %v", file, err)
						} else {
							if consFasta {
								fmt.Fprintf(out, "%60a\n", c)
							} else {
								fmt.Fprintf(out, "%q\n", c)
							}
							out.Close()
						}
					}
				}
			}()
		}
	}
	wait()
}

type familyReader struct {
	r       *gff.Reader
	last    *gff.Feature
	lastFam string
}

func (r *familyReader) readInto(v *[]*gff.Feature) error {
	*v = (*v)[:0]
	if r.last != nil {
		*v = append(*v, r.last)
	}
	for {
		f, err := r.r.Read()
		if err != nil {
			return err
		}
		gf := f.(*gff.Feature)
		fam := gf.FeatAttributes.Get("Family")
		isNext := fam != r.lastFam
		r.last = gf
		r.lastFam = fam
		if isNext && len(*v) != 0 {
			return nil
		}
		*v = append(*v, gf)
	}
}

func getReference(refName string) map[string]*linear.Seq {
	var f io.Reader
	f, err := os.Open(refName)
	if err != nil {
		log.Fatalf("error: %v", err)
	}
	defer f.(*os.File).Close()

	refStore := map[string]*linear.Seq{}
	if filepath.Ext(refName) == ".gz" {
		f, err = gzip.NewReader(f)
		if err != nil {
			log.Fatalf("failed to read reference: %v", err)
		}
	}
	r := fasta.NewReader(f, &linear.Seq{Annotation: seq.Annotation{Alpha: alphabet.DNA}})
	sc := seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		refStore[s.Name()] = s
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("failed to read reference: %v", err)
	}

	return refStore
}

var manager struct {
	limit chan struct{}
	wg    sync.WaitGroup
}

func acquire() {
	manager.wg.Add(1)
	manager.limit <- struct{}{}
}

func release() {
	<-manager.limit
	manager.wg.Done()
}

func wait() {
	manager.wg.Wait()
}

func consensus(in, aligner string) (*linear.QSeq, error) {
	var (
		m   *exec.Cmd
		err error
	)
	switch strings.ToLower(aligner) {
	case "muscle":
		m, err = muscle.Muscle{InFile: in, Quiet: true}.BuildCommand()
	case "mafft":
		m, err = mafft.Mafft{InFile: in, Auto: true, Quiet: true}.BuildCommand()
	default:
		log.Fatal("no valid aligner specified")
	}
	if err != nil {
		return nil, err
	}
	buf := &bytes.Buffer{}
	m.Stdout = buf
	err = m.Run()
	if err != nil {
		return nil, err
	}
	var (
		r  = fasta.NewReader(buf, &linear.Seq{Annotation: seq.Annotation{Alpha: alphabet.DNA}})
		ms = &multi.Multi{ColumnConsense: seq.DefaultQConsensus}
	)
	sc := seqio.NewScanner(r)
	for sc.Next() {
		ms.Add(sc.Seq())
	}
	return ms.Consensus(true), sc.Error()
}
