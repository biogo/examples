package main

import (
	"code.google.com/p/biogo/align/pals"
	"code.google.com/p/biogo/exp/alphabet"
	"code.google.com/p/biogo/exp/seq"
	"code.google.com/p/biogo/exp/seq/linear"
	"code.google.com/p/biogo/exp/seqio/fasta"
	"code.google.com/p/biogo/util"
	"crypto/md5"
	"fmt"
	"os"
	"path/filepath"
)

func packSequence(fileName string) (*pals.Packed, error) {
	_, name := filepath.Split(fileName)
	packer := pals.NewPacker(name)

	file, err := os.Open(fileName)
	if err == nil {
		md5hash, _ := util.Hash(md5.New(), file)
		logger.Printf("Reading %s: %s", fileName, fmt.Sprintf("%x", md5hash))

		template := &linear.Seq{Annotation: seq.Annotation{Alpha: alphabet.DNA}}
		seqFile := fasta.NewReader(file, template)

		f, p := logger.Flags(), logger.Prefix()
		if verbose {
			logger.SetFlags(0)
			logger.SetPrefix("")
			logger.Println("Sequence            \t    Length\t   Bin Range")
		}

		var seq seq.Sequence
		for {
			seq, err = seqFile.Read()
			if err != nil {
				break
			}
			s, err := packer.Pack(seq.(*linear.Seq))
			if err != nil {
				return nil, err
			}
			if verbose {
				logger.Println(s)
			}
		}
		if verbose {
			logger.SetFlags(f)
			logger.SetPrefix(p)
		}
	} else {
		logger.Fatalf("Error: %v.\n", err)
	}

	return packer.FinalisePack(), nil
}
