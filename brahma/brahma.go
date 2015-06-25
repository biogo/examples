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

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/util"
	"github.com/biogo/store/interval"
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
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	if *help || *sourceName == "" {
		flag.Usage()
		os.Exit(0)
	}

	if *targetName == "" {
		fmt.Fprintln(os.Stderr, "Reading PALS features from stdin.")
		target = gff.NewReader(os.Stdin)
	} else if tf, err := os.Open(*targetName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		fmt.Fprintf(os.Stderr, "Reading target features from `%s'.\n", *targetName)
		defer tf.Close()
		target = gff.NewReader(tf)
	}

	sf, err := os.Open(*sourceName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.\n", err)
		os.Exit(1)
	}
	fmt.Fprintf(os.Stderr, "Reading annotation features from `%s'.\n", *sourceName)
	defer sf.Close()
	source = gff.NewReader(sf)

	if *outName == "" {
		fmt.Fprintln(os.Stderr, "Writing annotation to stdout.")
		out = gff.NewWriter(os.Stdout, 60, false)
	} else if of, err := os.Create(*outName); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
	} else {
		defer of.Close()
		buf := bufio.NewWriter(of)
		defer buf.Flush()
		out = gff.NewWriter(buf, 60, true)
		fmt.Fprintf(os.Stderr, "Writing annotation to `%s'.\n", *outName)
	}
	out.Precision = 2

	ts := make(trees)

	var id uintptr
	for {
		f, err := source.Read()
		if err != nil {
			if err != io.EOF {
				fmt.Fprintf(os.Stderr, "Error: %v", err)
				os.Exit(1)
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
			fmt.Fprintf(os.Stderr, "missing repeat tag: file probably not an RM gff.\n")
			os.Exit(1)
		}
		repData.parse(ra)

		if t, ok := ts[gf.SeqName]; ok {
			err = t.Insert(repData, true)
		} else {
			t = &interval.IntTree{}
			err = t.Insert(repData, true)
			ts[gf.SeqName] = t
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, "insertion error: %v with repeat: %v\n", err, gf)
		}
	}
	for _, t := range ts {
		t.AdjustRanges()
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
				log.Fatalf("failed to read feature: %v", err)
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

		var consRemain int
		if consStart := rec.left; consStart != none {
			consEnd := rec.right
			consRemain = rec.remains
			repLen := consStart + consRemain
			if repLen <= 0 {
				repLen = util.MaxInt
			}

			if rec.genomic.Start() < target.FeatStart {
				consStart += target.FeatStart - rec.genomic.Start()
			}
			if rec.genomic.End() > target.FeatEnd {
				consEnd -= rec.genomic.End() - target.FeatEnd
			}

			leftMargin = float64(consStart) / float64(repLen)
			rightMargin = float64(consRemain) / float64(repLen)
		}

		mapStart := int(float64(start-target.FeatStart)*scale + 0.5)
		mapEnd := int(float64(end-target.FeatStart)*scale + 0.5)

		if target.FeatStrand == seq.Minus {
			mapStart, mapEnd = mapLen-mapEnd, mapLen-mapStart
		}

		if mapStart < 0 || mapStart > mapLen || mapEnd < 0 || mapEnd > mapLen {
			panic(fmt.Sprintf("brahma: failed to map: mapStart: %d, mapEnd: %d, mapLen: %d\n",
				mapStart, mapEnd, mapLen))
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
				mapping[mapEnd-1] &^= ('a' - 'A') // Uppercase has priority - truncation is indicated by fractional rep information.
			}
		}

		buf.WriteByte(' ')
		buf.WriteString(rec.name)

		if rec.left >= 0 {
			var (
				full    = float64(rec.right + consRemain)
				missing = float64(rec.left + consRemain)
			)
			if target.FeatStart > rec.genomic.Start() {
				missing += float64(target.FeatStart - rec.genomic.Start())
			}
			if rec.genomic.End() > target.FeatEnd {
				missing += float64(rec.genomic.End() - target.FeatEnd)
			}
			fmt.Fprintf(buf, "(%.0f%%)", ((full-missing)*100)/full)
		}
	}

	return buf.Bytes()
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

func (r *record) parse(a string) {
	fields := strings.Split(a, " ")

	r.name = fields[0]
	r.class = fields[1]
	if fields[2] != "." {
		r.left, _ = strconv.Atoi(fields[2])
	} else {
		r.left = none
	}
	if fields[3] != "." {
		r.right, _ = strconv.Atoi(fields[3])
	} else {
		r.right = none
	}
	if fields[4] != "." {
		r.remains, _ = strconv.Atoi(fields[4])
	} else {
		r.remains = none
	}
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
