source("~/scripts/tags.r")
source("~/scripts/peaks.r")

if(!interactive()) {
  args <- commandArgs(T)
  bamfile <- args[1]
  peakfile <- args[2]
}

pets <- readInBamFile(bamfile)
peaks <- read.narrow.peak.file(peakfile)

both.mapped <- !is.na(pets@pos) & !is.na(pets@mpos)

idx <- seq(1,sum(both.mapped),by=2)

first <- GRanges(seqnames=pets@rname[both.mapped][idx],
                 ranges=IRanges(pets@pos[both.mapped][idx],width=1),
                 strand='*')
second <- GRanges(seqnames=pets@rname[both.mapped][idx+1],
                 ranges=IRanges(pets@pos[both.mapped][idx+1],width=1),
                 strand='*')

of <- findOverlaps(first,peaks)
os <- findOverlaps(second,peaks)

counts <- matrix(0,nrow=length(peaks),ncol=length(peaks))

i <- intersect(queryHits(of),queryHits(os))

sf <- subjectHits(of)[queryHits(of) %in% i]
ss <- subjectHits(os)[queryHits(os) %in% i]

for( i in unique(sf) ) {
  w <- which(sf==i)
  n <- table(ss[w])
  
  counts[i,as.integer(names(n))] <- as.vector(n)
}

total <- rowSums(counts)
N <- sum(total)

lapply(1:(length(peaks)-1),function(i) {
  v <- c()
  for( j in (i+1):length(peaks)) {
    
    v <- c(v,i,j,)
  }
  
  matrix(v,nrow=length(v)/3,byrow=T)
})

#choose(total[2],counts[2,])[-2]*choose(2*N,total[-2]-counts[2,-2])/choose(2*N,total[-2])
#choose(total[51818],counts[51818,])[-51818]*choose(2*N-total[51818],total[-51818]-counts[51818,][-51818])/choose(2*N,total[-51818])
