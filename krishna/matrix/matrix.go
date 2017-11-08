// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// matrix runs a set of sequence segments (possibly chromosomes) through krisha
// performing self alignment on the diagonal and target/query alignment in the
// upper triangle.
package main

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"sync/atomic"
)

// Limit the number of parallel krishna jobs.
var (
	limit chan struct{}
	wg    sync.WaitGroup
	m     sync.Mutex
)

func acquire() {
	wg.Add(1)
	limit <- struct{}{}
}

func release(b *bytes.Buffer) {
	<-limit
	m.Lock()
	io.Copy(os.Stderr, b)
	m.Unlock()
	wg.Done()
}

var krishna string

func init() {
	var err error
	krishna, err = exec.LookPath("krishna")
	if err != nil {
		log.Fatalf("could not find krishna executable: %v", err)
	}
}

// Threadsafe counter.
var n int32

func done() int32 {
	return atomic.AddInt32(&n, 1)
}

func main() {
	kflags := flag.String("krishnaflags", "-tmp=/scratch -threads=2 -log", "Quoted set of flags to pass to krishna child processes.")
	threads := flag.Int("threads", 6, "Number of concurrent krishna instances to run.")
	flag.Parse()

	limit = make(chan struct{}, *threads)

	if len(flag.Args()) < 1 {
		log.Fatal("need targets")
	}
	files := flag.Args()
	t := (len(files)*len(files) + len(files)) / 2
	for _, f := range files {
		acquire()
		go runSelf(t, f, *kflags)
	}
	for i := range files[1:] {
		for j := range files[i : len(files)-1] {
			acquire()
			go runPair(t, files[i], files[j+i+1], *kflags)
		}
	}
	wg.Wait()
}

func runSelf(t int, target, kflags string) {
	b := &bytes.Buffer{}
	defer release(b)

	tbase := filepath.Base(target)
	if ext := filepath.Ext(tbase); len(ext) > 0 {
		tbase = tbase[:len(tbase)-len(ext)]
	}

	outfile := fmt.Sprintf("%s.gff", tbase)
	if _, err := os.Stat(outfile); err == nil {
		fmt.Fprintf(b, "file %q exists, skipping %d...\n", outfile, done())
		return
	}

	cmd := exec.Command(krishna, append(strings.Fields(kflags), "-target="+target, "-self", "-out="+outfile)...)
	cmd.Stderr = b
	err := cmd.Run()
	if err != nil {
		log.Printf("problem with %v self: %v\n", target, err)
	} else {
		b.Reset()
		fmt.Fprintf(b, "done %s, %d of %d\n", target, done(), t)
	}
}

func runPair(t int, target, query, kflags string) {
	b := &bytes.Buffer{}
	defer release(b)

	tbase := filepath.Base(target)
	if ext := filepath.Ext(tbase); len(ext) > 0 {
		tbase = tbase[:len(tbase)-len(ext)]
	}
	qbase := filepath.Base(query)
	if ext := filepath.Ext(qbase); len(ext) > 0 {
		qbase = qbase[:len(qbase)-len(ext)]
	}

	outfile := fmt.Sprintf("%s_%s.gff", tbase, qbase)
	if _, err := os.Stat(outfile); err == nil {
		fmt.Fprintf(b, "file %q exists, skipping %d...\n", outfile, done())
		return
	}

	cmd := exec.Command(krishna, append(strings.Fields(kflags), "-target="+target, "-query="+query, "-out="+outfile)...)
	cmd.Stderr = b
	err := cmd.Run()
	if err != nil {
		log.Printf("problem with %v and %v: %v\n", target, query, err)
	} else {
		b.Reset()
		fmt.Fprintf(b, "done %s x %s, %d of %d\n", target, query, done(), t)
	}
}
