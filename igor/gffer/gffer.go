// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// gffer converts the JSON output of igor to GFF.
package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"log"
	"os"

	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
)

type feat struct {
	Chr    string     `json:"C"`
	Start  int        `json:"S"`
	End    int        `json:"E"`
	Orient seq.Strand `json:"O"`
}

func main() {
	r := bufio.NewReader(os.Stdin)
	b := bufio.NewWriter(os.Stdout)
	defer b.Flush()
	w := gff.NewWriter(b, 60, false)

	ft := &gff.Feature{
		Source:         "igor",
		Feature:        "repeat",
		FeatAttributes: gff.Attributes{{Tag: "Family"}},
	}
	var v []*feat
	for fam := 0; ; fam++ {
		l, err := r.ReadBytes('\n')
		if err != nil {
			break
		}
		v = v[:0]
		err = json.Unmarshal(l, &v)
		if err != nil {
			log.Fatalf("error: %v", err)
		}
		for _, f := range v {
			ft.SeqName = f.Chr
			ft.FeatStart = f.Start
			ft.FeatEnd = f.End
			ft.FeatStrand = f.Orient
			ft.FeatFrame = gff.NoFrame
			ft.FeatAttributes[0].Value = fmt.Sprint(fam)
			_, err := w.Write(ft)
			if err != nil {
				log.Fatalf("error: %v", err)
			}
		}
	}
}
