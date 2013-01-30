// Copyright ©2011-2012 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package contig

import (
	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/seq"
	"code.google.com/p/biogo/seq/linear"

	"fmt"
	check "launchpad.net/gocheck"
	"testing"
)

func Test(t *testing.T) { check.TestingT(t) }

type S struct{}

var _ = check.Suite(&S{})

// Tests

type offsetSeq struct {
	seq    seq.Sequence
	offset int
}

type rep struct {
	val   string
	str   string
	fasta string
	rc    string
}

func (s *S) TestContig(c *check.C) {
	for i, t := range []struct {
		id      string
		length  int
		alpha   alphabet.Alphabet
		relaxed bool
		s       []offsetSeq
		rep     rep
		wrap    int
	}{
		{
			id:      "test",
			length:  20,
			alpha:   alphabet.DNA,
			relaxed: true,
			s: []offsetSeq{
				{
					linear.NewSeq("id1", alphabet.BytesToLetters([]byte("AGTC")), alphabet.DNA),
					2,
				},
				{
					linear.NewSeq("id2", alphabet.BytesToLetters([]byte("ACGT")), alphabet.DNA),
					15,
				},
			},
			rep: rep{
				val: `[0:n 2:"id1" AGTC 6:n 15:"id2" ACGT 19:n 20:<nil>]`, // If sequences overlap, this is misleading.
				str: `"test" nnAGTCnnnnnnnnnACGTn`,
				fasta: ">test\n" +
					"nnAGTCnnnnnnnnnACGTn",
				rc: `"test" nACGTnnnnnnnnnGACTnn`,
			},
		},
		{
			id:      "test",
			length:  20,
			alpha:   alphabet.DNA,
			relaxed: true,
			s: []offsetSeq{
				{
					linear.NewSeq("id1", alphabet.BytesToLetters([]byte("AGTC")), alphabet.DNA),
					2,
				},
				{
					linear.NewSeq("id2", alphabet.BytesToLetters([]byte("ACGT")), alphabet.DNA),
					15,
				},
			},
			rep: rep{
				val: `[0:n 2:"id1" AGTC 6:n 15:"id2" ACGT 19:n 20:<nil>]`, // If sequences overlap, this is misleading.
				str: `"test" nnAGTCnnnnnnnnnACGTn`,
				fasta: ">test\n" +
					"nnAGTCnnnn\n" +
					"nnnnnACGTn",
				rc: `"test" nACGTnnnnnnnnnGACTnn`,
			},
			wrap: 10,
		},
	} {
		con, err := New(t.id, t.length, t.alpha)
		c.Check(err, check.Equals, nil)
		con.Relaxed(t.relaxed)
		for _, os := range t.s {
			os.seq.SetOffset(os.offset)
			con.Insert(os.seq)
		}
		c.Check(fmt.Sprintf("%v", con), check.Equals, t.rep.val, check.Commentf("Test %d", i))
		c.Check(fmt.Sprintf("%s", con), check.Equals, t.rep.str, check.Commentf("Test %d", i))
		c.Check(fmt.Sprintf("%*a", t.wrap, con), check.Equals, t.rep.fasta, check.Commentf("Test %d", i))
		con.RevComp()
		c.Check(fmt.Sprintf("%s", con), check.Equals, t.rep.rc)
	}
}
