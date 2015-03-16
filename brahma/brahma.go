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
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"unsafe"

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

type Contig string

func (c Contig) Start() int             { return 0 }
func (c Contig) End() int               { return 0 }
func (c Contig) Len() int               { return 0 }
func (c Contig) Name() string           { return string(c) }
func (c Contig) Description() string    { return "Contig" }
func (c Contig) Location() feat.Feature { return nil }

type Location struct {
	Left, Right int
	Loc         feat.Feature
}

func (l Location) Start() int             { return l.Left }
func (l Location) End() int               { return l.Right }
func (l Location) Len() int               { return l.Right - l.Left }
func (l Location) Name() string           { return fmt.Sprintf("%s:[%d,%d)", l.Loc.Name(), l.Left, l.Right) }
func (l Location) Description() string    { return "Repeat" }
func (l Location) Location() feat.Feature { return l.Loc }

type RepeatRecord struct {
	Location

	Name, Class string
	Left, Right int
	Remain      int
}

func (rr *RepeatRecord) Overlap(b interval.IntRange) bool {
	return rr.Location.Right > b.Start && rr.Location.Left < b.End
}
func (rr *RepeatRecord) ID() uintptr { return uintptr(unsafe.Pointer(rr)) }
func (rr *RepeatRecord) Range() interval.IntRange {
	return interval.IntRange{rr.Location.Left, rr.Location.Right}
}

func (rr *RepeatRecord) Parse(a string) {
	fields := strings.Split(a, " ")

	rr.Name = fields[0]
	rr.Class = fields[1]
	if fields[2] != "." {
		rr.Left, _ = strconv.Atoi(fields[2])
	} else {
		rr.Left = -1
	}
	if fields[3] != "." {
		rr.Right, _ = strconv.Atoi(fields[3])
	} else {
		rr.Right = -1
	}
	if fields[4] != "." {
		rr.Remain, _ = strconv.Atoi(fields[4])
	} else {
		rr.Remain = -1
	}
}

type RepeatQuery struct {
	left, right int
	overlap     int
}

func (rq RepeatQuery) Overlap(b interval.IntRange) bool {
	return rq.right > b.Start+rq.overlap && rq.left < b.End-rq.overlap
}

type Match struct {
	Repeat  *RepeatRecord
	Overlap int
	Strand  seq.Strand
}

type Matches []Match

func (m Matches) Len() int {
	return len(m)
}
func (m Matches) Swap(i, j int) {
	m[i], m[j] = m[j], m[i]
}

type Overlap struct{ *Matches }

func (o Overlap) Less(i, j int) bool {
	return (*o.Matches)[i].Overlap < (*o.Matches)[j].Overlap
}
func (o Overlap) Pop() interface{} {
	*o.Matches = (*o.Matches)[:len(*o.Matches)-1]
	return nil
}
func (o Overlap) Push(x interface{}) {
	*o.Matches = append(*o.Matches, x.(Match))
}

type Start struct{ Matches }

func (s Start) Less(i, j int) bool {
	if s.Matches[i].Strand == 1 {
		return s.Matches[i].Repeat.Start() < s.Matches[j].Repeat.Start()
	}
	return s.Matches[i].Repeat.Start() > s.Matches[j].Repeat.Start()
}

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
	threads := flag.Int("threads", 2, "Number of threads to use.")
	bufferLen := flag.Int("buffer", 10, "Length of ouput buffer.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Parse()

	runtime.GOMAXPROCS(*threads)

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

	for {
		r, err := source.Read()
		if err != nil {
			if err != io.EOF {
				fmt.Fprintf(os.Stderr, "Error: %v", err)
				os.Exit(1)
			}
			break
		}
		repeat := r.(*gff.Feature)
		repData := &RepeatRecord{
			Location: Location{
				Left:  repeat.FeatStart,
				Right: repeat.FeatEnd,
				Loc:   Contig(repeat.SeqName),
			},
		}
		ra := repeat.FeatAttributes.Get("Repeat")
		if ra == "" {
			fmt.Fprintf(os.Stderr, "Missing repeat tag: file probably not an RM gff.\n")
			os.Exit(1)
		}
		repData.Parse(ra)

		if t, ok := ts[repeat.SeqName]; ok {
			err = t.Insert(repData, true)
		} else {
			t = &interval.IntTree{}
			err = t.Insert(repData, true)
			ts[repeat.SeqName] = t
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, "Insertion error: %v with repeat: %v\n", err, repeat)
		}
	}
	for _, t := range ts {
		t.AdjustRanges()
	}

	process := make(chan *gff.Feature)
	buffer := make(chan *gff.Feature, *bufferLen)
	processWg := &sync.WaitGroup{}
	outputWg := &sync.WaitGroup{}

	if *threads < 2 {
		*threads = 2
	}
	for i := 0; i < *threads-1; i++ {
		processWg.Add(1)
		go processServer(ts, process, buffer, processWg)
	}

	//output server
	outputWg.Add(1)
	go func() {
		defer outputWg.Done()
		for feature := range buffer {
			out.Write(feature)
		}
	}()

	for {
		feature, err := target.Read()
		if err == nil {
			process <- feature.(*gff.Feature)
		} else {
			close(process)
			break
		}
	}

	processWg.Wait()
	close(buffer)
	outputWg.Wait()
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

var blankTemplate = `"` + strings.Repeat("-", mapLen)

func processServer(index trees, queue, output chan *gff.Feature, wg *sync.WaitGroup) {
	defer wg.Done()
	const tag = "Annot"
	var (
		blank = []byte(blankTemplate)

		buffer  = make([]byte, 0, annotationLength)
		mapping = buffer[1 : mapLen+1]
		annots  = make(Matches, 0, maxAnnotations+1)
		o       = Overlap{&annots}
		overlap int
	)

	for f := range queue {
		overlap = int(float64(f.Len()) * minOverlap)
		annots = annots[:0] // Obviates heap initialisation.
		buffer = buffer[:len(blank)]
		copy(buffer, blank)

		t, ok := index[f.SeqName]
		if ok {
			t.DoMatching(func(hit interval.IntInterface) (done bool) {
				r := hit.Range()
				heap.Push(o, Match{
					Repeat:  hit.(*RepeatRecord),
					Overlap: min(r.End, f.FeatEnd) - max(r.Start, f.FeatStart),
					Strand:  f.FeatStrand,
				})
				if len(annots) > maxAnnotations {
					heap.Pop(o)
				}
				return
			}, RepeatQuery{f.FeatStart, f.FeatEnd, overlap})
		}

		if len(annots) > 1 {
			sort.Sort(Start{annots})
		}
		if len(annots) > 0 {
			buffer = makeAnnot(f, annots, mapping, bytes.NewBuffer(buffer))
		}

		buffer = append(buffer, '"')
		f.FeatAttributes = append(f.FeatAttributes, gff.Attribute{
			Tag:   tag,
			Value: string(buffer),
		})

		output <- f
	}
}

func makeAnnot(f *gff.Feature, m Matches, mapping []byte, buf *bytes.Buffer) []byte {
	var leftMargin, rightMargin float64
	scale := float64(mapLen) / float64(f.Len())
	for i, ann := range m {
		var (
			rep      = ann.Repeat
			location = rep.Location
			start    = max(location.Start(), f.FeatStart)
			end      = min(location.End(), f.FeatEnd)
		)

		var consRemain = 0
		if consStart := rep.Left; consStart != -1 {
			consEnd := rep.Right
			consRemain = rep.Remain
			repLen := consStart + consRemain
			if repLen <= 0 {
				repLen = util.MaxInt
			}

			if location.Start() < f.FeatStart {
				consStart += f.FeatStart - location.Start()
			}
			if location.End() > f.FeatEnd {
				consEnd -= location.End() - f.FeatEnd
			}

			leftMargin = float64(consStart) / float64(repLen)
			rightMargin = float64(consRemain) / float64(repLen)
		}

		mapStart := int(float64(start-f.FeatStart)*scale + 0.5)
		mapEnd := int(float64(end-f.FeatStart)*scale + 0.5)

		if f.FeatStrand == -1 {
			mapStart, mapEnd = mapLen-mapEnd, mapLen-mapStart
		}

		if mapStart < 0 || mapStart > mapLen || mapEnd < 0 || mapEnd > mapLen {
			panic(fmt.Sprintf("brahma: failed to map: mapStart: %d, mapEnd: %d, mapLen: %d\n",
				mapStart, mapEnd, mapLen))
		}

		if mapStart < mapEnd {
			cLower := 'a' + byte(i)
			cUpper := 'A' + byte(i)

			if leftMargin <= maxMargin && rep.Left != -1 {
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
				if rightMargin <= maxMargin && rep.Left != -1 {
					mapping[mapEnd-1] = cUpper
				} else {
					mapping[mapEnd-1] = cLower
				}
			} else if rightMargin <= maxMargin && rep.Left != -1 {
				mapping[mapEnd-1] &^= ('a' - 'A') // Uppercase has priority - truncation is indicated by fractional rep information.
			}
		}

		buf.WriteByte(' ')
		buf.WriteString(rep.Name)

		if rep.Left >= 0 {
			var (
				full    = float64(rep.Right + consRemain)
				missing = float64(rep.Left + consRemain)
			)
			if f.FeatStart > location.Start() {
				missing += float64(f.FeatStart - location.Start())
			}
			if location.End() > f.FeatEnd {
				missing += float64(location.End() - f.FeatEnd)
			}
			fmt.Fprintf(buf, "(%.0f%%)", ((full-missing)*100)/full)
		}
	}

	return buf.Bytes()
}
