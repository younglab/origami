library(compiler,quietly = !interactive())
invisible(enableJIT(3)) ## enable JIT compilation

if( !interactive() ) {
  args <- commandArgs(T)
  peakcounts <- if( !is.na(args[1])) args[1] else "peak-counts.txt"
  intcounts <- if( !is.na(args[2])) args[2] else "int-counts.txt"
  outfile <- if( !is.na(args[3])) args[3] else "results.csv"
  modelfile <- if( !is.na(args[4])) args[4] else "model-data.Rdata"
  iterations <- if( !is.na(args[5])) as.integer(args[5]) else 10000
  burnin <- if( !is.na(args[6])) as.integer(args[6]) else 100
  prune <- if( !is.na(args[7])) as.integer(args[7]) else 5
  minimodel <- if( !is.na(args[8])) args[8]=="yes" else T
  usedistance <- if( !is.na(args[9])) args[9]=="yes" else T
  usedf <- if( !is.na(args[10])) as.integer(args[10]) else 0
  useglm <- if( !is.na(args[11])) args[11] == "yes" else F
  
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

b <- i1 != i2


p <- intcounts[i1 != i2,]

#totint <- tapply(rep(p$V7,2),factor(c(as.character(i1[b]),as.character(i2[b])),levels=levels(i1)),mean)
totint <- table(factor(c(as.character(i1[b]),as.character(i2[b])),levels=levels(i1)))

inttable <- table(factor(c(as.character(i1[i1!=i2]),as.character(i2[i1!=i2])),levels=levels(i1)),rep(p$V7,2))

if( !interactive() ){
  
  cat("Running hypergeometric test...\n")
  hyperg <- estimate.hypergeometric.pvalue(p,depth)

  cat("Running two-component Bayesian mixture model...\n")

    #cat("Model 1\n")
  gbayes.m1 <- estimate.global.bayesian.mixture(p,depth,inttable,burnin=burnin,N=iterations+burnin,pruning=prune,show.progress=T,mini.model=minimodel,
                                                with.distance.weight=usedistance,usedf=usedf,useglm=useglm)
  gbayesp1 <- extract.global.bayesian.mixture.prob(gbayes.m1)
  
  #cat("Model 2\n")
  #gbayes.m2 <- estimate.global.bayesian.weighed.depth.mixture(p,depth,burnin=burnin,N=iterations+burnin,pruning=prune,show.progress=T)
  #gbayesp2 <- extract.global.bayesian.mixture.prob(gbayes.m2)


  
  m <- cbind(p,hyperg,gbayesp1)#,gbayesp2)
  
  colnames(m) <- c("chromosome1","start1","end1","chromosome2","start2","end2","PET Count","Hypergeometric p-value",
                   "Bayes mixture 1")#,"Bayes mixture 2")

  write.csv(m,file=outfile,row.names=F,quote=F)
  save(hyperg,gbayes.m1,file=modelfile)
  
  #save(hyperg,gbayes.m1,gbayes.m2,file=modelfile)
}