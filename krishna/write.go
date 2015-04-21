package main

import (
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"os"
	"sync"
	"unsafe"

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
	tf, err := os.Create(fmt.Sprintf("%s-%s.traps.le.gz", outFile, d))
	if err != nil {
		return err
	}
	gz := gzip.NewWriter(tf)
	// TODO(kortschak): Write int size to file so we are arch independent.
	err = binary.Write(gz, binary.LittleEndian, unsafeTraps(traps))
	if err != nil {
		return err
	}
	err = gz.Close()
	if err != nil {
		return err
	}
	return tf.Close()
}

func init() {
	switch unsafe.Sizeof(int(0)) {
	case unsafe.Sizeof(int64(0)), unsafe.Sizeof(int32(0)):
	default:
		panic("int type unknown size")
	}
}

func unsafeTraps(traps []filter.Trapezoid) interface{} {
	switch unsafe.Sizeof(int(0)) {
	case unsafe.Sizeof(int64(0)):
		type trapezoid64 struct{ Top, Bottom, Left, Right int64 }
		return *(*[]trapezoid64)(unsafe.Pointer(&traps))
	case unsafe.Sizeof(int32(0)):
		type trapezoid32 struct{ Top, Bottom, Left, Right int32 }
		return *(*[]trapezoid32)(unsafe.Pointer(&traps))
	default:
		panic("int type unknown size")
	}
}
