library(GenomicRanges)

estimate.hypergeometric.pvalue <- function(ints,depth) {
  d <- GRanges(seqnames=as.character(depth$V1),ranges=IRanges(depth$V2,depth$V3),strand='*')
  g1 <- GRanges(seqnames=as.character(ints$V1),ranges=IRanges(ints$V2,ints$V3),strand='*')
  g2 <- GRanges(seqnames=as.character(ints$V4),ranges=IRanges(ints$V5,ints$V6),strand='*')
  
  m1 <- match(g1,d)
  m2 <- match(g2,d)
  pets <- ints[,7]
  dv <- depth[,4]
  pval <- c()
  pval <- lchoose(dv[m1],pets)+lchoose(2*sum(dv)-dv[m1],dv[m2]-pets) - lchoose(2*sum(dv),dv[m2])
  
  exp(pval)
}

#choose(total[2],counts[2,])[-2]*choose(2*N,total[-2]-counts[2,-2])/choose(2*N,total[-2])
#choose(total[51818],counts[51818,])[-51818]*choose(2*N-total[51818],total[-51818]-counts[51818,][-51818])/choose(2*N,total[-51818])
