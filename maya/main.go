// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// maya reads a collection of motif features from a BED file and finds motifs
// in this collection that fall within regions specified in a second BED file.
// The mean and variance of motif locations of motifs found found is reported.
package main

import (
	"code.google.com/p/biogo.interval"
	"code.google.com/p/biogo/io/featio/bed"
	"flag"
	"fmt"
	"math"
	"os"
	"unsafe"
)

type trees map[string]*interval.IntTree

type Motif struct {
	Start, End int
	Contig     string
}

func (m *Motif) Overlap(b interval.IntRange) bool {
	return m.End > b.Start && m.Start < b.End
}
func (m *Motif) ID() uintptr              { return uintptr(unsafe.Pointer(m)) }
func (m *Motif) Range() interval.IntRange { return interval.IntRange{m.Start, m.End} }
func (m *Motif) String() string           { return fmt.Sprintf("%s\t%d\t%d", m.Contig, m.Start, m.End) }

type Region struct {
	Start, End int
	Contig     string
}

func (r *Region) Overlap(b interval.IntRange) bool {
	return r.Start <= b.Start && b.End <= r.End
}
func (r *Region) ID() uintptr              { return uintptr(unsafe.Pointer(r)) }
func (r *Region) Range() interval.IntRange { return interval.IntRange{r.Start, r.End} }
func (r *Region) String() string           { return fmt.Sprintf("%s\t%d\t%d", r.Contig, r.Start, r.End) }

func main() {
	motifName := flag.String("motif", "", "Filename for motif file.")
	regionName := flag.String("region", "", "Filename for region file.")
	verbose := flag.Bool("verbose", false, "Print details of identified motifs to stderr.")
	headerLine := flag.Bool("header", false, "Print a header line.")
	help := flag.Bool("help", false, "Print this usage message.")

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: %s -motif <motif file> -region <region file>\n", os.Args[0])
		flag.PrintDefaults()
	}

	flag.Parse()

	if *help || *regionName == "" || *motifName == "" {
		flag.Usage()
		os.Exit(1)
	}

	// Open files
	motifFile, err := os.Open(*motifName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(0)
	}
	defer motifFile.Close()
	motif := bed.NewReader(motifFile, 3)
	fmt.Fprintf(os.Stderr, "Reading motif features from `%s'.\n", *motifName)

	regionFile, err := os.Open(*regionName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(0)
	}
	defer regionFile.Close()
	region := bed.NewReader(regionFile, 3)
	fmt.Fprintf(os.Stderr, "Reading region features from `%s'.\n", *regionName)

	// Read in motif features and build interval tree to search
	ts := make(trees)

	for line := 1; ; line++ {
		motifLine, err := motif.Read()
		if err != nil {
			break
		}

		motif := &Region{
			Start:  motifLine.Start(),
			End:    motifLine.End(),
			Contig: fmt.Sprint(motifLine.Location()),
		}
		if t, ok := ts[motif.Contig]; ok {
			err = t.Insert(motif, true)
		} else {
			t = &interval.IntTree{}
			err = t.Insert(motif, true)
			ts[motif.Contig] = t
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, "Insertion error: %v with motif: %v\n", err, motif)
		}

	}

	// Read in region features and search for motifs within region
	// Calculate median motif location, sample standard deviation of locations
	// and mean distance of motif from midpoint of region for motifs contained
	// within region. Report these and n of motifs within region.
	if *headerLine {
		fmt.Println("Chromosome\tStart\tEnd\tn-hits\tMeanHitPos\tStddevHitPos\tMeanMidDistance")
	}
	for line := 1; ; line++ {
		regionLine, err := region.Read()
		if err != nil {
			break
		} else {
			if regionLine.Start() > regionLine.End() {
				fmt.Fprintf(os.Stderr, "Line: %d: Feature has end < start: %v\n", line, regionLine)
				continue
			}
			region := &Region{
				Start:  regionLine.Start(),
				End:    regionLine.End(),
				Contig: fmt.Sprint(regionLine.Location()),
			}
			regionMidPoint := float64(region.Start+region.End) / 2
			if *verbose {
				fmt.Fprintln(os.Stderr, region)
			}
			sumOfDiffs, sumOfSquares, mean, oldmean, n := 0., 0., 0., 0., 0.

			if t, ok := ts[region.Contig]; ok {
				t.DoMatching(func(m interval.IntInterface) (done bool) {
					r := m.Range()
					mid := float64(r.Start+r.End) / 2
					if *verbose {
						fmt.Fprintf(os.Stderr, "\t%s\n", m)
					}

					// The Method of Provisional Means	
					n++
					mean = oldmean + (mid-oldmean)/n
					sumOfSquares += (mid - oldmean) * (mid - mean)
					oldmean = mean

					sumOfDiffs += math.Abs(mid - regionMidPoint)

					return
				}, region)
			}
			fmt.Printf("%s\t%d\t%d\t%0.f\t%0.f\t%f\t%f\n",
				region.Contig, region.Start, region.End,
				n, mean, math.Sqrt(sumOfSquares)/(n-1), sumOfDiffs/n)
		}
	}
}
