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

if(length(args)<6) {
  stop("args: binding file, int file, N PETs, P(self), P(random)")
}

binding <- read.in.binding.sites(args[1])
ints <- read.in.interactions(args[2])
npets <- as.integer(args[3])
pself <- as.numeric(args[4])
prandom <- as.numeric(args[5])
fout <- args[6]
pint <- 1-pself-prandom
print(pint)

#if(sum(c(pself,prandom,pint)!=1)) {
#  stop("probabilities need to add up to 1")
#}

m <- do.call(rbind,replicate(npets,simulate.pets(),simplify = F))

z <- table(m[,1],m[,2])

bs <- paste(seqnames(binding),start(binding),end(binding),sep='\t')


#####

#write.table(z,file=fout,quote=F,row.names=F,col.names=F,sep='\t')

df <- data.frame(bs,sapply(1:nrow(z),function(i) sum(c(z[i,],z[,i],-z[i,i]))))

write.table(df,file=paste(fout,".a",sep=''),quote=F,row.names=F,col.names=F,sep='\t')

df <- data.frame(s1=NULL,s2=NULL,cn=NULL)

for( i in 1:(nrow(z)-1) ) {
  for( j in (i+1):nrow(z)) {
    df <- rbind(df,data.frame(s1=bs[i],s2=bs[j],cn=z[i,j]))
  }
}

write.table(df,file=paste(fout,".b",sep=''),quote=F,row.names=F,col.names=F,sep='\t')


#for(i in 1:nrow(z)) {
  
#}

