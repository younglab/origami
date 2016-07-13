library(Rsamtools)

args <- commandArgs(T)

if(length(args)<2) {
  stop("make-chip-profile.r <bam file> <output file>")
}

bfile <- args[1]
ofile <- args[2]

param <- ScanBamParam(what=c("rname","pos"),flag=scanBamFlag(isUnmappedQuery=F))
d <- scanBam(bfile,param=param)[[1]]

ntotal <- length(d$pos)

cn <- lapply(split(d$pos,d$rname),function(v) {
  if(length(v)==0) return(NULL)
  r <- range(v)
  breaks <- seq(r[1],r[2]+50,by=50)
  
  list(r[1],hist(v,breaks=breaks,plot=F)$counts)
})

f <- file(ofile,"w")

cat("track type=wiggle_0 name=\"ChIA-PET\" description=\"ChIA-PET density\"",file = f,sep='\n')

dummy <- mapply(function(l,n) {
  if(is.null(l[[1]])) return(NULL)
  cat(paste("fixedStep chrom=",n," start=",format(l[[1]],scientific = F)," step=50",sep=''),file=f,sep='\n')
  writeLines(format(l[[2]],scientific = F),con=f,sep='\n')
},cn,as.list(names(cn)),SIMPLIFY=F)

close(f)
