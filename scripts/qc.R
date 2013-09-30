library(qrqc)
library(ggplot2)
library(gridExtra)

pdf("results/qc.pdf")

the_args <- commandArgs(TRUE)
# the_args <- c("test.fastq", "test.fastq", "test.fastq", "test.fastq", "test.fastq")

makeGraph<-function(fn) {
  s<-strsplit(fn, "=")
  label<-s[[1]][[1]]
  fn<-s[[1]][[2]]
  qualPlot(readSeqFile(fn, quality='sanger', hash=FALSE)) + 
      opts(title = label) +
      opts(plot.title = theme_text(size = 10)) +
      opts(axis.title.x=theme_blank(), axis.title.y=theme_blank())
}

plots<-lapply(the_args, makeGraph)

do.call(grid.arrange, c(plots, list(ncol=4)))

dev.off()

