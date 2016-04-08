args <- commandArgs(T)

peakfile <- args[1]
outfile <- args[2]
threshold <- as.integer(args[3])

peaks <- read.table(peakfile,sep='\t')

counts <- log10(peaks$V4)

if(any(is.infinite(counts))) {
  b <- is.infinite(counts)
  counts <- counts[!b]
  warning("\tSome peaks were removed with having 0 PET counts, this generally happens when using a set of peak calls not derived from the ChIA-PET data")
}

zscore <- (counts - mean(counts))/sd(counts)

m <- peaks[zscore<threshold,1:3]

write.table(m,file=outfile,sep='\t',quote=F,row.names=F,col.names=F)