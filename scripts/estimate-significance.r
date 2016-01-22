source("~/dsday/origami/scripts/hypergeometric-test.r")
#source("~/dsday/origami/scripts/per-sample-bayesian.r")
source("~/dsday/origami/scripts/estimate-global-bayesian-mixture.r")

peakcounts <- "peak-counts.txt"
intcounts <- "int-counts.txt"
outfile <- "results.csv"

if( !interactive() ) {
  args <- commandArgs(T)
  if( !is.na(args[1])) peakcounts <- args[1]
  if( !is.na(args[2])) intcounts <- args[2]
  if( !is.na(args[3])) outfile <- args[3]
}

convert.to.factor <- function(pos,idx,l=NULL) {
  d <- pos[,idx]
  s <- paste(d[,1],d[,2],d[,3],sep='_')
  
  r <- if(is.null(l)) factor(s) else factor(s,levels=levels(l)) 
  r
}


depth <- read.table(peakcounts,sep='\t')
intcounts <- read.table(intcounts,sep='\t')

f <- convert.to.factor(depth,1:3)

i1 <- convert.to.factor(intcounts,1:3,f)
i2 <- convert.to.factor(intcounts,4:6,f)

p <- intcounts[i1 != i2,]

if( !interactive() ){
  
  cat("Running hypergeometric test...\n")
  hyperg <- estimate.hypergeometric.pvalue(p,depth)
  #bayesps <- estimate.per.sample.bayesian.probability(p,depth)
  
  cat("Running two-component Bayesian mixture model...\n")
  gbayes.m <- estimate.global.bayesian.mixture(p,depth,show.progress=T)
  gbayesp <- extract.global.bayesian.prob(gbayes.m)
  
  #m <- cbind(p,hyperg,bayesps,gbayesp)
  m <- cbind(p,hyperg,gbayesp)
  
  #colnames(m) <- c("chromosome1","start1","end1","chromosome2","start2","end2","PET Count","Hypergeometric p-value","Bayes posterior probability 1","Bayes global mixture posterior probability")
  colnames(m) <- c("chromosome1","start1","end1","chromosome2","start2","end2","PET Count","Hypergeometric p-value","Bayes global mixture posterior probability")
  
  write.csv(m,file=outfile,row.names=F,quote=F)
}