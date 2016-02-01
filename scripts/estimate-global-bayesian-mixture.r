library(utils)
library(matrixStats)

estimate.global.bayesian.mixture <- function(ints,depth,N=1100,burnin=100,pruning=NULL,with.distance.weight=F,no.depth=F,multiply=F,show.progress=F) {
  S <- nrow(ints)
  
  d <- GRanges(seqnames=as.character(depth$V1),ranges=IRanges(depth$V2,depth$V3),strand='*')
  g1 <- GRanges(seqnames=as.character(ints$V1),ranges=IRanges(ints$V2,ints$V3),strand='*')
  g2 <- GRanges(seqnames=as.character(ints$V4),ranges=IRanges(ints$V5,ints$V6),strand='*')
  
  m1 <- match(g1,d)
  m2 <- match(g2,d)
  counts <- ints[,7]
  d <- depth[,4]
  
  #mdepth <- mean(apply(cbind(m1,m2),1,function(v) sum(d[v])))
  if(!multiply) {
    sdepth <- rowSums(cbind(d[m1],d[m2]))
    msdepth <- mean(sdepth)
  } else {
    sdepth <- rowProds(cbind(d[m1],d[m2]))
    msdepth <- mean(sdepth)
  }
  
  print(c(msdepth,range(sdepth)))
  
  l <- lapply(1:S,function(i) list(z=rep(NA_integer_,N),p1=c(.5,rep(NA_real_,N))))
  pp <- rep(.5,S)
  
  lambda0 <- c(1,rep(NA_real_,N))
  lambda1 <- c(5,rep(NA_real_,N))

  
  totcounts <- sum(counts)

  
  if(show.progress) pb <- txtProgressBar()
    
  
  #f <- function(p,h1,h0) {
  #  g1 <- p*dpois(counts,h1)
  #  g2 <- (1-p) * dpois(counts,h0)
  #  
  #  g1/rowSums(cbind(g1,g2)) 
  #}
  
  for( i in 1:N ) {
    #vp <- f(pp,lambda1[i],lambda0[i])
    g1 <- pp*dpois(counts,lambda1[i])
    g2 <- (1-pp) * dpois(counts,lambda0[i])
    
    vp <- g1/rowSums(cbind(g1,g2)) 
    if(any(is.na(vp))) vp[is.na(vp)] <- 0 ### need  more intelligent way to handle this
    vz <- rbinom(S,1,vp)

    pp <- if( no.depth ) rbeta(S,1+vz,1+(1-vz)) else rbeta(S,sdepth/msdepth+vz,msdepth/msdepth+(1-vz))
    
    l <- mapply(function(lx,z,p){
      lx$z[i] <- z
      lx$p1[i+1] <- p
      lx
    },l,as.list(vz),as.list(pp),SIMPLIFY=F)
    
    l0 <- 0
    r <- sum(counts[vz==0])
    n <- sum(vz==0)
    
    while(l0 <= 0 ) l0 <- rgamma(1,r,n)

    l1 <- l0
    #r <- totcounts - r #sum(counts[vz==1])

    l1x <- rgamma(1,totcounts-r,S-n)#sum(vz==1))
    l1 <- max(l1,l1x)
    
    lambda0[i+1] <- l0
    lambda1[i+1] <- l1

    print(c(l0,l1))
    if(show.progress) setTxtProgressBar(pb,i/N)
  }
  
  if(show.progress) close(pb)
  ret <- list(s=l,l0=lambda0,l1=lambda1)
  if(!is.null(burnin) && is.integer(burnin) && burnin > 1) {
    orig <- ret
    idx <- -(1:burnin)
    ret <- lapply(ret,function(v) v[idx])
    ret$orig <- orig
  }
  ret
}

extract.global.bayesian.prob <- function(model) {
  sapply(model$s,function(v) sum(v$z)/length(v$z))
}
