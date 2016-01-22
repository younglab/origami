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
  
  mdepth <- mean(apply(cbind(m1,m2),1,function(v) sum(d[v])))
  
  l <- lapply(1:nrow(ints),function(i) list(sdepth=d[m1[i]]+d[m2[i]],z=c(),p1=c(.5)))
  
  lambda0 <- c(1)
  lambda1 <- c(5)
  hsd0 <- hsd1 <- sd(counts)
  p1 <- c(.5)
  #a0 <- c(1)
  #b0 <- c(1)
  #a1 <- c(1)
  #b1 <- c(1)
  
  if(show.progress) pb <- txtProgressBar()
    
  
  f <- function(p,h1,h0) p*dpois(counts,h1)/rowSums(cbind(p*dpois(counts,h1),(1-p)*dpois(counts,h0)))
  
  for( i in 1:N ) {
#    print(c(lambda1[i],lambda0[i]))
    vp <- f(sapply(l,function(v) v$p1[i]),lambda1[i],lambda0[i])
    if(any(is.na(vp))) vp[is.na(vp)] <- 0 ### need  more intelligent way to handle this
    vz <- rbinom(S,1,vp)
    #print(sum(vz))
    l <- mapply(function(lx,z){
      p <- if( no.depth ) rbeta(1,1+z,1+(1-z)) else rbeta(1,lx$sdepth+z,mdepth+(1-z))
      lx$z[i] <- z
      lx$p1[i+1] <- p
      lx
    },l,as.list(vz),SIMPLIFY=F)
    
    l0 <- 0
    r <- sum(counts[vz==0])
    #print(r)
   # print((a[i]*lambda0[i]+sum(vz==0)*r)/(a[i]+sum(vz)))
    while(l0 <= 0 ) l0 <- rgamma(1,r,sum(vz==0))
    #a0[i+1] <- a0[i]+r
    #b0[i+1] <- b0[i]+sum(vz==0)
#    print(l0)
    l1 <- l0
    r <- sum(counts[vz==1])
    #print(r)
    
    #print(c(a1[i]+r,b1[i]+sum(vz==1)))
    while(l1 <= l0 ) {
      l1 <- rgamma(1,r,sum(vz==1))

      #print(c(l0,l1))
    }
    #a1[i+1] <- a1[i]+r
    #b1[i+1] <- b1[i]+sum(vz==1)
    
    
    lambda0[i+1] <- l0
    lambda1[i+1] <- l1
    #print('---')
    if(show.progress) setTxtProgressBar(pb,i/N)
  }
  
  if(show.progress) close(pb)
  list(s=l,l0=lambda0,l1=lambda1)#,a0=a0,b0=b0,a1=a1,b1=b1)
}

extract.global.bayesian.prob <- function(model) {
  sapply(model$s,function(v) sum(v$z)/length(v$z))
}
