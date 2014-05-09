// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package turner

import (
	"code.google.com/p/biogo.store/interval"
	"code.google.com/p/biogo.store/interval/landscape"
	"code.google.com/p/biogo/align/pals"
)

// Lambda is a collection of lambda functions for a single sequence position.
type Lambda []int

func (l Lambda) at(k int) int {
	if k < len(l) {
		return l[k]
	}
	return 0
}

type ivs []interval.IntRange

func (s ivs) Len() int                     { return len(s) }
func (s ivs) Less(i, j int) bool           { return s[i].Start < s[j].Start }
func (s ivs) Swap(i, j int)                { s[i], s[j] = s[j], s[i] }
func (s ivs) Item(i int) interval.IntRange { return s[i] }

// Landscape holds persistence landscape calculation results produced by Paint.
type Landscape struct {
	Chromosome string
	Start, End int

	// Note contains any additional information.
	Note string

	// Lambdas holds the set of lambda functions for each position within the pile.
	Lambdas []Lambda

	// MaxK is the maximum depth of the landscape.
	MaxK int

	// Features is a representation of edge features in the landscape.
	Features [][]float64
}

// condMax returns the max of a and b if cond == true, otherwise a.
func condMax(cond bool, a, b int) int {
	if cond && b > a {
		return b
	}
	return a
}

// Paint caluclates a persistence landscape data set for rendering.
//
// An example renderer is provided in the following R snippet which takes the
// JSON-marshaled serialisation of a Landscape:
//
/*
	require(jsonlite)
	library(grid)
	library(gridExtra)
	library(ggplot2)
	library(reshape)
	library(spatstat)

	fill <- function (l, maxK) {
		matrix(unlist(lapply(l, function(x) { x[1:maxK] })), ncol=length(l), nrow=maxK)
	}

	turner <- function(file, doblur=F, sigma=1, bleed=TRUE, thresh=0.1, adjust=0.1) {
		a <- fromJSON(file)

		lambdas <- melt(fill(a$Lambdas, a$MaxK))
		names(lambdas) <-c ("k", "pos", "lambda")
		p1 <- ggplot(lambdas, aes(x=pos, y=k, fill=lambda)) +
			xlab(NULL) +
			geom_raster() +
			scale_fill_gradient(high = "white") +
			theme(line = element_blank(), panel.background = element_blank())
		g1 <- ggplot_gtable(ggplot_build(p1))


		features <- fill(a$Features, a$MaxK)
		if(doblur) {
			features <- melt(as.matrix(blur(
				im(features, 1:ncol(features), 1:nrow(features)),
				sigma=sigma, bleed=bleed)))
		} else {
			features <- melt(features)
		}
		names(features)<-c("k", "pos", "delta")

		features.c <- rowSums(xtabs(delta ~ pos+k, features[ which(features$delta >= 0), ]))
		features.c <- data.frame(cbind(as.numeric(names(features.c)), features.c))
		names(features.c) <- c("pos", "crest")
		features.t <- rowSums(xtabs(delta ~ pos+k, features[ which(features$delta <= 0), ]))
		features.t <- data.frame(cbind(as.numeric(names(features.t)), features.t))
		names(features.t) <- c("pos", "trough")
		both <- melt(merge(features.c, features.t, by="pos", all=TRUE), id="pos")
		both[is.na(both)] <- 0
		names(both) <- c("pos", "feature", "score")

		features <- features[ which(abs(features$delta) > thresh), ]
		p2 <- ggplot(features, aes(x=pos, y=k, fill=delta)) +
			xlab(NULL) +
			geom_raster() +
			scale_fill_gradient2(low = "red", mid="white") +
			theme(line = element_blank(), panel.background = element_blank())
		g2 <- ggplot_gtable(ggplot_build(p2))

		p3 <- ggplot(data=both, aes(x=pos, y=score, color=feature)) +
			xlab(NULL) +
			geom_line()
		g3 <- ggplot_gtable(ggplot_build(p3))

		g2$widths <- g1$widths
		g3$widths <- g1$widths

		grid.arrange(g1, g2, g3, nrow = 3,
			main=sprintf("turner of %s:%d-%d\n%s", a$Chromosome, a$Start, a$End, a$Note))
	}

	turner("landscape-chr11-22306794-22312360.json")
*/
func Paint(p *pals.Pile, all bool) Landscape {
	data := make(ivs, 0, len(p.Images))
	for _, im := range p.Images {
		data = append(data, interval.IntRange{Start: im.Start(), End: im.End()})
	}

	ls := Landscape{
		Start: p.From,
		End:   p.To,

		// We describe ends as well, so we need length+1 for lambdas.
		Lambdas:  make([]Lambda, p.Len()+1),
		Features: make([][]float64, p.Len()-1),
	}
	landscape.Describe(data, func(pos int, l []int) {
		ls.Lambdas[pos-p.From] = append(Lambda(nil), l...)
		if len(l) > ls.MaxK {
			ls.MaxK = len(l)
		}
	})
	if all {
		for i, pos := range ls.Lambdas {
			ls.Lambdas[i] = append(pos, make([]int, ls.MaxK-len(pos))...)
		}
	}

	for i, pos := range ls.Lambdas[1 : len(ls.Lambdas)-1] {
		var pf []float64
		for k := 0; k < condMax(all, len(pos), ls.MaxK); k++ {
			var delta float64
			h := pos.at(k)
			switch {
			// Trough detection.
			case h == 0 && ls.Lambdas[i].at(k)+ls.Lambdas[i+2].at(k) > 0,
				(i == 0 || i == len(ls.Lambdas)-3) && h > 0,
				(h <= ls.Lambdas[i].at(k) && h < ls.Lambdas[i+2].at(k)),
				(h < ls.Lambdas[i].at(k) && h <= ls.Lambdas[i+2].at(k)):
				delta = -1
				if h != 0 {
					delta /= float64(h)
				}
			// Peak detection.
			case (h >= ls.Lambdas[i].at(k) && h > ls.Lambdas[i+2].at(k)),
				(h > ls.Lambdas[i].at(k) && h >= ls.Lambdas[i+2].at(k)):
				delta = 1
			}
			pf = append(pf, delta)
		}
		ls.Features[i] = pf
	}

	return ls
}
