library(compiler,quietly = !interactive())
invisible(enableJIT(3)) ## enable JIT compilation

if( !interactive() ) {
  args <- commandArgs(T)
  peakcounts <- if( !is.na(args[1])) args[1] else "peak-counts.txt"
  intcounts <- if( !is.na(args[2])) args[2] else "int-counts.txt"
  outfile <- if( !is.na(args[3])) args[3] else "results.csv"
  modelfile <- if( !is.na(args[4])) args[4] else "model-data.Rdata"
  
  args <- commandArgs()
  
  f <- sub("--file=","",args[grep("--file=",args)])
  
  dbase <- dirname(f)
}

source(paste(dbase,"hypergeometric-test.r",sep='/'))
source(paste(dbase,"estimate-global-bayesian-mixture.r",sep='/'))



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

  cat("Running two-component Bayesian mixture model...\n")
  gbayes.m <- estimate.global.bayesian.mixture(p,depth,show.progress=T)
  gbayesp <- extract.global.bayesian.prob(gbayes.m)
  
  gbayesnd.m <- estimate.global.bayesian.no.depth.mixture(p,depth,show.progress=T)
  gbayesndp <- extract.global.bayesian.prob(gbayesnd.m)
  
  m <- cbind(p,hyperg,gbayesp,gbayesndp)
  
  colnames(m) <- c("chromosome1","start1","end1","chromosome2","start2","end2","PET Count","Hypergeometric p-value",
                   "Bayes global mixture posterior probability","Bayes No Depth Mixture Posterior Probability")

  write.csv(m,file=outfile,row.names=F,quote=F)
  save(hyperg,gbayes.m,gbayesnd.m,file=modelfile)
}