library(GenomicRanges)

read.in.binding.sites <- function(fname) {
  temp <- read.table(fname)
  
  GRanges(seqnames=temp$V1,ranges=IRanges(temp$V2,temp$V3),strand='*')
}

read.in.interactions <- function(fname) {
  temp <- read.table(fname)
  
  temp[,3] <- temp[,3]/sum(temp[,3])
  
  temp
}

simulate.pets <- function() {
  p <- runif(1,0,1)
  
  if(p<prandom) { ## random ligation
    i <- sort(sample.int(length(binding),size=2,replace=T))
    r <- c(i[1],i[2])
  } else if(prandom < p && p < (prandom+pself)) {
    idx <- sample.int(length(binding),size=1)
    r <- c(idx,idx)
  } else {
    idx <- sample.int(nrow(ints),size=1,prob=ints[,3])
    r <- c(ints[idx,1],ints[idx,2])
  }
  
  r
}

####

if( !interactive() ) args <- commandArgs(T)

if(length(args)<5) {
  stop("args: binding file, int file, N PETs, P(self), P(random)")
}

binding <- read.in.binding.sites(args[1])
ints <- read.in.interactions(args[2])
npets <- as.integer(args[3])
pself <- as.numeric(args[4])
prandom <- as.numeric(args[5])
pint <- 1-pself-prandom
print(pint)

#if(sum(c(pself,prandom,pint)!=1)) {
#  stop("probabilities need to add up to 1")
#}

m <- do.call(rbind,replicate(npets,simulate.pets(),simplify = F))

z <- table(m[,1],m[,2])


#####