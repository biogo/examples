package main

import (
	"code.google.com/p/biogo/align/pals"
	"code.google.com/p/biogo/align/pals/dp"
	"sync"
)

var wlock = &sync.Mutex{}

func WriteDPHits(w *pals.Writer, target, query *pals.Packed, hits []dp.DPHit, comp bool) (n int, err error) {
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
