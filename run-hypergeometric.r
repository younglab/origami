if(!interactive()) args <- commandArgs(T)

m <- read.table(args[1],sep='\t',comment.char="")

b <- duplicated(m$V7)
idx <- sort(c(w,w-1))

pl <- unique(msub$V14)
msub <- m[idx,]
counts <- matrix(0,nrow=length(pl),ncol=length(pl))
rownames(counts) <- colnames(counts) <- pl

pc <- as.character(msub$V14)

midx <- match(pc,pl)

i <- 1
while(i < nrow(msub)) {
  j <- i + 1
  
  counts[midx[i],midx[j]] <- counts[midx[i],midx[j]] + 1
  i <- i+2
}