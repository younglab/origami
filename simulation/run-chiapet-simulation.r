library(GenomicRanges)

read.in.binding.sites <- function(fname) {
  temp <- read.table(fname,sep='\t')
  
  GRanges(seqnames=temp$V1,ranges=IRanges(temp$V2,temp$V3),strand='*')
}

read.in.interactions <- function(fname) {
  temp <- read.table(fname,sep='\t')
  
  temp
}

args <- commandArgs(T)

binding <- read.in.binding.sites(args[1])
ints <- read.in.interactions(args[2])


