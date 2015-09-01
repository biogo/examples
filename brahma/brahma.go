// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// brahma performs annotation of GFF intervals produced by PALS/PILER, taking
// annotation information from a GFF file generated from RepeatMasker output.
package main

import (
	"bufio"
	"bytes"
	"container/heap"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/store/interval"
	"github.com/biogo/store/step"
)

const (
	maxAnnotations   = 8
	annotationLength = 256
	mapLen           = 20
	maxMargin        = 1. / mapLen
)

var (
	minOverlap float64
)

type trees map[string]*interval.IntTree

func main() {
	var (
		target *gff.Reader
		source *gff.Reader
		out    *gff.Writer
		err    error
	)

	targetName := flag.String("target", "", "Filename for input to be annotated. Defaults to stdin.")
	sourceName := flag.String("source", "", "Filename for source annotation.")
	outName := flag.String("out", "", "Filename for output. Defaults to stdout.")
	flag.Float64Var(&minOverlap, "overlap", 0.05, "Overlap between features.")
	covRep := flag.String("covrep", "", "Filename for repeat type coverage report.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	if *help || *sourceName == "" {
		flag.Usage()
		os.Exit(0)
	}

	if *targetName == "" {
		fmt.Fprintln(os.Stderr, "reading PALS features from stdin.")
		target = gff.NewReader(os.Stdin)
	} else if tf, err := os.Open(*targetName); err != nil {
		log.Fatalf("could not open %q: %v", *targetName, err)
	} else {
		fmt.Fprintf(os.Stderr, "reading target features from %q.\n", *targetName)
		defer tf.Close()
		target = gff.NewReader(tf)
	}

	sf, err := os.Open(*sourceName)
	if err != nil {
		log.Fatalf("could not open %q: %v", *sourceName, err)
	}
	fmt.Fprintf(os.Stderr, "reading annotation features from %q.\n", *sourceName)
	defer sf.Close()
	source = gff.NewReader(sf)

	if *outName == "" {
		fmt.Fprintln(os.Stderr, "writing annotation to stdout.")
		out = gff.NewWriter(os.Stdout, 60, false)
	} else if of, err := os.Create(*outName); err != nil {
		log.Fatalf("could not create %q: %v", err)
	} else {
		defer of.Close()
		buf := bufio.NewWriter(of)
		defer buf.Flush()
		out = gff.NewWriter(buf, 60, true)
		fmt.Fprintf(os.Stderr, "writing annotation to %q.\n", *outName)
	}
	out.Precision = 2

	ts := make(trees)

	var id uintptr
	for {
		f, err := source.Read()
		if err != nil {
			if err != io.EOF {
				log.Fatalf("failed to read source feature: %v", err)
			}
			break
		}

		gf := f.(*gff.Feature)
		repData := &record{
			id: id,
			genomic: repeat{
				left:  gf.FeatStart,
				right: gf.FeatEnd,
				loc:   contig(gf.SeqName),
			},
		}
		id++

		ra := gf.FeatAttributes.Get("Repeat")
		if ra == "" {
			log.Fatal("missing repeat tag: file probably not an RM gff.")
		}
		err = repData.parse(ra)
		if err != nil {
			log.Fatalf("failed to parse repeat tag: %v\n", err)
		}

		if t, ok := ts[gf.SeqName]; ok {
			err = t.Insert(repData, true)
		} else {
			t = &interval.IntTree{}
			err = t.Insert(repData, true)
			ts[gf.SeqName] = t
		}
		if err != nil {
			log.Fatalf("insertion error: %v with repeat: %v\n", err, gf)
		}
	}
	for _, t := range ts {
		t.AdjustRanges()
	}

	var coverage map[string][2]*step.Vector
	if *covRep != "" {
		coverage = make(map[string][2]*step.Vector)
	}

	const tag = "Annot"
	var (
		blank = `"` + strings.Repeat("-", mapLen)

		buffer  = make([]byte, 0, annotationLength)
		mapping = buffer[1 : mapLen+1]
		annots  = make(matches, 0, maxAnnotations+1)
		best    = byOverlap{&annots}
		overlap int
	)
	for {
		rf, err := target.Read()
		if err != nil {
			if err != io.EOF {
				log.Fatalf("failed to read target feature: %v", err)
			}
			break
		}
		f := rf.(*gff.Feature)

		overlap = int(float64(f.Len()) * minOverlap)
		annots = annots[:0] // Obviates heap initialisation.
		buffer = buffer[:len(blank)]
		copy(buffer, blank)

		t, ok := ts[f.SeqName]
		if ok {
			t.DoMatching(func(hit interval.IntInterface) (done bool) {
				r := hit.Range()
				heap.Push(best, match{
					record:  hit.(*record),
					overlap: min(r.End, f.FeatEnd) - max(r.Start, f.FeatStart),
					strand:  f.FeatStrand,
				})
				if len(annots) > maxAnnotations {
					// byOverlap is a min heap for overlap,
					// so pop removes the lowest overlap.
					heap.Pop(best)
				}
				return
			}, query{f.FeatStart, f.FeatEnd, overlap})
		}

		if len(annots) > 1 {
			sort.Sort(byStart{annots})
		}
		if *covRep != "" {
			for _, a := range annots {
				if a.record.left == none {
					continue
				}
				v, ok := coverage[a.record.name]
				if !ok {
					for i := range v {
						v[i], err = step.New(0, a.record.right+a.record.remains, stepBool(false))
						if err != nil {
							panic(err)
						}
						v[i].Relaxed = true // This should not be required, but RepeatMasker.
					}
					coverage[a.record.name] = v
				}

				// krishna coverage.
				left := a.record.left + max(0, f.FeatStart-a.record.genomic.Start())
				right := a.record.right + min(0, f.FeatEnd-a.record.genomic.End())
				if right < left { // This craziness is... because RepeatMasker.
					continue
				}
				v[1].SetRange(left, right, stepBool(true))
			}
		}

		if len(annots) > 0 {
			buffer = makeAnnot(f, annots, mapping, bytes.NewBuffer(buffer))
		}

		buffer = append(buffer, '"')
		f.FeatAttributes = append(f.FeatAttributes, gff.Attribute{
			Tag:   tag,
			Value: string(buffer),
		})

		out.Write(f)
	}

	if *covRep != "" {
		// RepeatMasker coverage for repeat types seen by krishna.
		for _, t := range ts {
			t.Do(func(iv interval.IntInterface) (done bool) {
				rec := iv.(*record)
				if rec.left == none || rec.right < rec.left { // RepeatMasker...
					return
				}
				if v, ok := coverage[rec.name]; ok {
					v[0].SetRange(rec.left, rec.right, stepBool(true))
				}
				return
			})
		}
		err = writeCoverage(*covRep, coverage)
		if err != nil {
			log.Fatalf("failed to write coverage report: %v", err)
		}
	}
}

// stepBool is a bool type satisfying the step.Equaler interface.
type stepBool bool

// Equal returns whether b equals e. Equal assumes the underlying type of e is a stepBool.
func (b stepBool) Equal(e step.Equaler) bool {
	return b == e.(stepBool)
}

// makeAnoot return an annotation map
func makeAnnot(target *gff.Feature, m matches, mapping []byte, buf *bytes.Buffer) []byte {
	var leftMargin, rightMargin float64
	scale := mapLen / float64(target.Len())
	for i, annotation := range m {
		var (
			rec   = annotation.record
			start = max(rec.genomic.Start(), target.FeatStart)
			end   = min(rec.genomic.End(), target.FeatEnd)
		)

		// Handle defined length repeat margins.
		var consRemain int
		if rec.left != none {
			consStart := rec.left + max(0, target.FeatStart-rec.genomic.Start())
			consEnd := rec.right + min(0, target.FeatEnd-rec.genomic.End())
			consRemain = rec.remains
			length := consEnd + consRemain

			leftMargin = float64(consStart) / float64(length)
			rightMargin = float64(consEnd) / float64(length)
		}

		mapStart := int(float64(start-target.FeatStart)*scale + 0.5)
		mapEnd := int(float64(end-target.FeatStart)*scale + 0.5)

		if target.FeatStrand == seq.Minus {
			mapStart, mapEnd = mapLen-mapEnd, mapLen-mapStart
		}

		if mapStart < 0 || mapStart > mapLen || mapEnd < 0 || mapEnd > mapLen {
			log.Fatalf("failed to map: mapStart: %d, mapEnd: %d, mapLen: %d feature: %+v annotation: %+v\n",
				mapStart, mapEnd, mapLen, target, rec)
		}

		if mapStart < mapEnd {
			cLower := 'a' + byte(i)
			cUpper := 'A' + byte(i)

			if leftMargin <= maxMargin && rec.left != none {
				mapping[mapStart] = cUpper
			} else {
				mapping[mapStart] = cLower
			}

			if mapEnd-mapStart > 1 {
				seg := mapping[mapStart+1 : mapEnd-1]
				for p := range seg {
					seg[p] = cLower
				}
			}

			if mapEnd-1 != mapStart {
				if rightMargin <= maxMargin && rec.left != none {
					mapping[mapEnd-1] = cUpper
				} else {
					mapping[mapEnd-1] = cLower
				}
			} else if rightMargin <= maxMargin && rec.left != none {
				// Uppercase has priority - truncation is
				// indicated by fractional rep information.
				mapping[mapEnd-1] &^= ('a' - 'A')
			}
		}

		buf.WriteByte(' ')
		buf.WriteString(rec.name)

		if rec.left != none {
			fmt.Fprintf(buf, "(%.0f%%|%.0f%%)",
				// Overlap with masked element.
				float64(annotation.overlap)/float64(rec.genomic.Len())*100,
				// Overlap with complete element.
				float64(annotation.overlap)/float64(rec.right+rec.remains)*100,
			)
		}
	}

	return buf.Bytes()
}

func writeCoverage(file string, coverage map[string][2]*step.Vector) error {
	if len(coverage) == 0 {
		return nil
	}
	f, err := os.Create(file)
	if err != nil {
		return err
	}
	defer f.Close()

	const mapLen = 40
	cov := [2][]byte{make([]byte, mapLen), make([]byte, mapLen)}

	tw := tabwriter.NewWriter(f, 0, 0, 1, ' ', 0)
	fmt.Fprintln(tw, "repeat\tnamed coverage\tde novo coverage")
	defer tw.Flush()

	names := make([]string, 0, len(coverage))
	for n := range coverage {
		names = append(names, n)
	}
	sort.Strings(names)
	for _, n := range names {
		pair := coverage[n]
		scale := mapLen / float64(max(pair[0].Len(), pair[1].Len()))
		for i, v := range pair {
			for j := range cov[0] {
				cov[i][j] = '-'
			}
			v.Do(func(start, end int, e step.Equaler) {
				if e.(stepBool) {
					for j := int(float64(start) * scale); j < int(float64(end)*scale); j++ {
						cov[i][j] = '#'
					}
				}
			})
		}
		fmt.Fprintf(tw, "%s\t%s\t%s\n", n, cov[0], cov[1])
	}

	return nil
}

// contig is a sequence contig with repeats mapped to it.
type contig string

func (c contig) Start() int             { return 0 }
func (c contig) End() int               { return 0 }
func (c contig) Len() int               { return 0 }
func (c contig) Name() string           { return string(c) }
func (c contig) Description() string    { return "contig" }
func (c contig) Location() feat.Feature { return nil }

// repeat is a repeat-matching interval.
type repeat struct {
	left, right int
	loc         feat.Feature
}

func (r repeat) Start() int             { return r.left }
func (r repeat) End() int               { return r.right }
func (r repeat) Len() int               { return r.right - r.left }
func (r repeat) Name() string           { return fmt.Sprintf("%s:[%d,%d)", r.loc.Name(), r.left, r.right) }
func (r repeat) Description() string    { return "repeat" }
func (r repeat) Location() feat.Feature { return r.loc }

// record is a masked repeat record.
type record struct {
	id uintptr

	// genomic is the genomic region matched
	// to the the repeat identified below.
	genomic repeat

	// name and class that the repeat type
	// and class defined by the masker.
	name, class string

	// left and right are the left and right
	// position of the record alignment in
	// consensus-relative coordinates.
	left, right int
	// remains is the distance from right to
	// the end of the consensus sequence.
	remains int
}

func (r *record) Overlap(b interval.IntRange) bool {
	return r.genomic.End() > b.Start && r.genomic.Start() < b.End
}
func (r *record) ID() uintptr { return r.id }
func (r *record) Range() interval.IntRange {
	return interval.IntRange{r.genomic.Start(), r.genomic.End()}
}

const none = -1

func (r *record) parse(a string) error {
	fields := strings.Split(a, " ")

	r.name = fields[0]
	r.class = fields[1]
	var err error
	if fields[2] != "." {
		r.left, err = strconv.Atoi(fields[2])
		if err != nil {
			return err
		}
		r.left-- // Convert to 0-based indexing.
	} else {
		r.left = none
	}
	if fields[3] != "." {
		r.right, err = strconv.Atoi(fields[3])
		if err != nil {
			return err
		}
	} else {
		r.right = none
	}
	if fields[4] != "." {
		r.remains, err = strconv.Atoi(fields[4])
		if err != nil {
			return err
		}
	} else {
		r.remains = none
	}

	return nil
}

// query is an interval query allowing for an overlap threshold.
type query struct {
	left, right int
	overlap     int
}

func (q query) Overlap(b interval.IntRange) bool {
	return q.right > b.Start+q.overlap && q.left < b.End-q.overlap
}

// match is a repeat match to a target interval.
type match struct {
	// record is the matching repeat record.
	record *record

	// overlap is the overlap between the
	// target and the record.
	overlap int

	// strand is the strand of the target.
	strand seq.Strand
}

type matches []match

func (m matches) Len() int {
	return len(m)
}
func (m matches) Swap(i, j int) {
	m[i], m[j] = m[j], m[i]
}

// byOverlap is a min heap of matches by overlap.
type byOverlap struct{ *matches }

func (o byOverlap) Less(i, j int) bool {
	return (*o.matches)[i].overlap < (*o.matches)[j].overlap
}
func (o byOverlap) Pop() interface{} {
	*o.matches = (*o.matches)[:len(*o.matches)-1]
	return nil
}
func (o byOverlap) Push(x interface{}) {
	*o.matches = append(*o.matches, x.(match))
}

// byStart implements the sort.Interface, sorting by start position.
type byStart struct{ matches }

func (s byStart) Less(i, j int) bool {
	if s.matches[i].strand == seq.Plus {
		return s.matches[i].record.genomic.Start() < s.matches[j].record.genomic.Start()
	}
	return s.matches[i].record.genomic.Start() > s.matches[j].record.genomic.Start()
}

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
