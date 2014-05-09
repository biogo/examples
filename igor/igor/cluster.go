// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package igor

import (
	"code.google.com/p/biogo.examples/igor/turner"

	"code.google.com/p/biogo.store/interval"
	"code.google.com/p/biogo.store/step"
	"code.google.com/p/biogo/align/pals"

	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"sort"
	"sync"
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

func overlap(a, b interval.IntRange) int {
	return max(0, max(a.End-b.Start, b.End-a.Start))
}

type pileInterval struct {
	p  *pals.Pile
	id uintptr
}

func (pi pileInterval) ID() uintptr { return pi.id }
func (pi pileInterval) Overlap(b interval.IntRange) bool {
	return pi.p.Start() < b.End && pi.p.End() > b.Start
}
func (pi pileInterval) Range() interval.IntRange {
	return interval.IntRange{pi.p.Start(), pi.p.End()}
}

// ClusterConfig specifies Cluster behaviour.
type ClusterConfig struct {
	// BandWidth specifies the maximum fractional distance between
	// endpoints of images being clustered into sub-piles. See
	// turner.Cluster for details.
	BandWidth float64

	// RequiredCover specifies the target coverage fraction for
	// each input pile. If RequiredCover is greater than 1, all
	// all sub-piles are retained depending on the values of
	// KeepOverlaps and OverlapThresh.
	RequiredCover float64

	// KeepOverlaps allows overlapping piles to be retained
	// during clustering. If false, contained features overlapping
	// by more than OverlapThresh fraction of the smaller pile are
	// are discarded.
	KeepOverlaps  bool
	OverlapThresh float64

	// LandscapeDir specifies the path to store persistence
	// landscape data. No data is stored if empty.
	LandscapeDir string

	// Threads specifies the number of independent clustering
	// instances to run in parallel. If zero, only single threaded
	// operation is performed.
	Threads int
}

// threadManager limits the number of runable clustering threads and
// handles logging control.
type threadManager struct {
	loggers chan *logger
	m       sync.Mutex
	wg      sync.WaitGroup
}

func (m *threadManager) acquire() *logger {
	m.wg.Add(1)
	return <-m.loggers
}

func (m *threadManager) release(l *logger) {
	m.m.Lock()
	io.Copy(os.Stderr, l.buf)
	m.m.Unlock()
	m.loggers <- l
	m.wg.Done()
}

func (m *threadManager) wait() {
	m.wg.Wait()
}

type logger struct {
	buf *bytes.Buffer
	log *log.Logger
}

func (l *logger) printf(format string, args ...interface{}) {
	l.log.Printf(format, args...)
}

// Cluster performs sub-pile clustering according to the config provided.
// The number of sub-piles and a collection of piles broken into sub-piles is
// returned.
func Cluster(piles []*pals.Pile, cfg ClusterConfig) (int, [][]*pals.Pile) {
	m := threadManager{
		loggers: make(chan *logger, max(cfg.Threads, 1)),
	}
	for len(m.loggers) < cap(m.loggers) {
		buf := &bytes.Buffer{}
		m.loggers <- &logger{
			buf: buf,
			log: log.New(buf, "", log.LstdFlags),
		}
	}

	clust := make([][]*pals.Pile, len(piles))
	// skipLock protect writes/reads to p.Loc which is abused as a flag to
	// allow Group to know which piles to ignore in the grouping phase.
	var skipLock sync.Mutex
	for i, p := range piles {
		i, p := i, p

		l := m.acquire()
		go func() {
			defer m.release(l)

			skipLock.Lock()
			loc := p.Loc
			skipLock.Unlock()
			if loc == nil {
				return
			}

			tc := turner.Cluster(p, cfg.BandWidth)
			logLine := fmt.Sprintf("%s:%d-%d; turner found %d clusters from %d members",
				loc.Name(), p.Start(), p.End(), len(tc), len(p.Images),
			)
			l.printf(logLine)

			sv, err := step.New(p.Start(), p.End(), step.Int(0))
			if err != nil {
				panic(err)
			}
			sort.Sort(turner.ByDepth(tc))
			var (
				t        interval.IntTree
				accepted int
			)
			for j, c := range tc {
				if !cfg.KeepOverlaps {
					pi := pileInterval{c, uintptr(j)}
					for _, iv := range t.Get(pi) {
						r := iv.Range()
						if within(cfg.OverlapThresh, overlap(r, pi.Range()), min(pi.p.Len(), r.End-r.Start)) {
							c = nil
							pile := iv.(pileInterval).p
							skipLock.Lock()
							pile.Loc = nil
							for _, im := range pile.Images {
								im.Location().(*pals.Pile).Loc = nil
							}
							skipLock.Unlock()
						}
					}
					if c == nil {
						tc[j] = nil
						continue
					}
					t.Insert(pi, false)
				}

				accepted++
				l.printf("cluster %d: %d members with %d volume spanning %d-%d",
					j, len(c.Images), turner.Volume(c), c.Start(), c.End(),
				)
				sv.ApplyRange(c.Start(), c.End(), func(e step.Equaler) step.Equaler {
					return e.(step.Int) + step.Int(len(c.Images))
				})
				var cov int
				sv.Do(func(start, end int, e step.Equaler) {
					if e.(step.Int) > 1 {
						cov += end - start
					}
				})
				if !within(cfg.RequiredCover, cov, p.Len()) {
					skipLock.Lock()
					for _, dc := range tc[j+1:] {
						dc.Loc = nil
						for _, im := range dc.Images {
							im.Location().(*pals.Pile).Loc = nil
						}
					}
					skipLock.Unlock()
					tc = tc[:j+1]
					break
				}
			}

			clust[i] = tc

			if cfg.LandscapeDir != "" && accepted > 5 {
				l.printf("writing turner landscape for %s:%d-%d", loc.Name(), p.From, p.To)
				ls := turner.Paint(p, false)
				ls.Chromosome = loc.Name()
				ls.Note = logLine
				err := os.Mkdir(filepath.Join(cfg.LandscapeDir, ls.Chromosome), 0755)
				if err != nil && !os.IsExist(err) {
					l.printf("failed to create subdirectory for : %q error: %v", ls.Chromosome, err)
					return
				}
				path := filepath.Join(cfg.LandscapeDir, ls.Chromosome, fmt.Sprintf("landscape-%s-%d-%d.json", ls.Chromosome, ls.Start, ls.End))
				f, err := os.Create(path)
				if err != nil {
					l.printf("failed to create landscape file: %q error: %v", path, err)
					return
				}
				defer f.Close()
				err = json.NewEncoder(f).Encode(ls)
				if err != nil {
					l.printf("failed to marshal turner painting: %v", err)
					return
				}
			}
		}()
	}
	m.wait()

	var n int
	for _, c := range clust {
		n += len(c)
	}
	return n, clust
}
