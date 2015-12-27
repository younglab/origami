source("~/dsday/origami/scripts/hypergeometric-test.r")
source("~/dsday/origami/scripts/per-sample-bayesian.r")


convert.to.factor <- function(pos,idx,l=NULL) {
  d <- pos[,idx]
  s <- paste(d[,1],d[,2],d[,3],sep='_')
  
  r <- if(is.null(l)) factor(s) else factor(s,levels=levels(l)) 
  r
}


depth <- read.table("peak-counts.txt",sep='\t')
intcounts <- read.table("int-counts.txt",sep='\t')

f <- convert.to.factor(depth,1:3)

i1 <- convert.to.factor(intcounts,1:3,f)
i2 <- convert.to.factor(intcounts,4:6,f)

p <- intcounts[i1 != i2,]

if( !interactive() ){
  hyperg <- estimate.hypergeometric.pvalue(p,depth)
  bayesps <- estimate.per.sample.bayesian.probability(p,depth)
  
  write.csv(cbind(p,hyperg,bayesps),file='results.csv')
}