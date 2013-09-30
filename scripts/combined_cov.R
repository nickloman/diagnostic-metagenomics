#!/usr/bin/env Rscript

library(GenomeGraphs)
library(GenomicFeatures)
library(ShortRead)
library(HilbertVis)
library(ggplot2)
library(gridExtra)

pdf("combined_depth.pdf")

args<-commandArgs(TRUE)

makeGraph<-function(fn) {
	aln<-readGappedAlignments(fn)
	cov<-coverage(aln)
	covnum<-as.numeric(cov$scaffold00001)
	covframe<-data.frame(pos=seq(0, length(covnum), length.out=5000), cov=shrinkVector(as.vector(covnum), 5000, mode="mean"))
	plot1<-ggplot(covframe, aes(pos, cov), log="y") + geom_point(alpha=1/10) + stat_smooth() + scale_y_continuous(limits=c(0,60))
}

plots<-lapply(args, makeGraph)
do.call(grid.arrange, c(plots, list(nrow=1, ncol=length(args))))

dev.off()

#ggsave(file = args[2], plot1, width=4.8, height=3.6)
#ggsave(file = paste(args[2], sep="", ".png"), plot1, width=4.8, height=3.6)

