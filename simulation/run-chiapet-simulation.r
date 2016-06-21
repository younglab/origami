library(GenomicRanges)

read.in.binding.sites <- function(fname) {
  temp <- read.table(fname,sep='\t')
  
  GRanges(seqnames=temp$V1,ranges=IRanges(temp$V2,temp$V3),strand='*')
}

read.in.interactions <- function(fname) {
  temp <- read.table(fname,sep='\t')
  
  temp
}

simulate.pets <- function() {
  p <- runif(1,0,1)
  
  if(p>.8) {
    i <- sort(sample.int(length(binding),size=2,replace=T))
    r <- c(i[1],i[2])
  } else {
    idx <- sample.int(nrow(ints),size=1)
    r <- c(ints[idx,1],ints[idx,2])
  }
  
  r
}

####

if( !interactive() ) args <- commandArgs(T)

binding <- read.in.binding.sites(args[1])
ints <- read.in.interactions(args[2])
npets <- as.integer(args[3])

m <- do.call(rbind,replicate(npets,simulate.pets(),simplify = F))

z <- table(m[,1],m[,2])


#####