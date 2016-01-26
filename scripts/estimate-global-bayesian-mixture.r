library(utils)

estimate.global.bayesian.mixture <- function(ints,depth,N=1000,burnin=NULL,pruning=NULL,with.distance.weight=F,no.depth=F,show.progress=F) {
  S <- nrow(ints)
  
  d <- GRanges(seqnames=as.character(depth$V1),ranges=IRanges(depth$V2,depth$V3),strand='*')
  g1 <- GRanges(seqnames=as.character(ints$V1),ranges=IRanges(ints$V2,ints$V3),strand='*')
  g2 <- GRanges(seqnames=as.character(ints$V4),ranges=IRanges(ints$V5,ints$V6),strand='*')
  
  m1 <- match(g1,d)
  m2 <- match(g2,d)
  counts <- ints[,7]
  d <- depth[,4]
  
  #mdepth <- mean(apply(cbind(m1,m2),1,function(v) sum(d[v])))
  sdepth <- rowSums(cbind(d[m1],d[m2]))
  msdepth <- mean(sdepth)
  
  l <- lapply(1:S,function(i) list(z=c(),p1=c(.5)))
  pp <- rep(.5,nrow(ints))
  
  lambda0 <- c(1)
  lambda1 <- c(5)

  
  if(show.progress) pb <- txtProgressBar()
    
  
  f <- function(p,h1,h0) {
    g1 <- p*dpois(counts,h1)
    g2 <- (1-p) * dpois(counts,h0)
    
    g1/rowSums(cbind(g1,g2))
  }
  
  for( i in 1:N ) {
    vp <- f(pp,lambda1[i],lambda0[i])
    if(any(is.na(vp))) vp[is.na(vp)] <- 0 ### need  more intelligent way to handle this
    vz <- rbinom(S,1,vp)

    pp <- if( no.depth ) rbeta(S,1+vz,1+(1-vz)) else rbeta(S,sdepth+vz,msdepth+(1-vz))
    
    l <- mapply(function(lx,z,p){
      lx$z[i] <- z
      lx$p1[i+1] <- p
      lx
    },l,as.list(vz),as.list(pp),SIMPLIFY=F)
    
    l0 <- 0
    r <- sum(counts[vz==0])
    
    while(l0 <= 0 ) l0 <- rgamma(1,r,sum(vz==0))

    l1 <- l0
    r <- sum(counts[vz==1])

    l1x <- rgamma(1,r,sum(vz==1))
    l1 <- max(l1,l1x)
    
    lambda0[i+1] <- l0
    lambda1[i+1] <- l1

    if(show.progress) setTxtProgressBar(pb,i/N)
  }
  
  if(show.progress) close(pb)
  list(s=l,l0=lambda0,l1=lambda1)#,a0=a0,b0=b0,a1=a1,b1=b1)
}

extract.global.bayesian.prob <- function(model) {
  sapply(model$s,function(v) sum(v$z)/length(v$z))
}
