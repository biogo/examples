package main

import (
	"encoding/gob"
	"fmt"
	"os"
	"sync"

	"github.com/biogo/biogo/align/pals"
	"github.com/biogo/biogo/align/pals/dp"
	"github.com/biogo/biogo/align/pals/filter"
)

var wlock = &sync.Mutex{}

func WriteDPHits(w *pals.Writer, target, query *pals.Packed, hits []dp.Hit, comp bool) (n int, err error) {
	wlock.Lock()
	defer wlock.Unlock()

	for _, hit := range hits {
		pair, err := pals.NewPair(target, query, hit, comp)
		if err != nil {
			return n, err
		} else {
			ln, err := w.Write(pair)
			n += ln
			if err != nil {
				return n, err
			}
		}
	}

	return
}

func WriteTraps(comp bool, traps filter.Trapezoids) error {
	var d string
	if comp {
		d = "rev"
	} else {
		d = "fwd"
	}
	tf, err := os.Create(fmt.Sprintf("%s-%s.gob", outFile, d))
	if err != nil {
		return err
	}
	enc := gob.NewEncoder(tf)
	err = enc.Encode(traps)
	if err != nil {
		return err
	}
	return tf.Close()
}
