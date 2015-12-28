estimate.global.bayesian.mixture <- function(ints,depth,N=1000) {
  S <- nrow(ints)
  
  l <- lapply(1:nrow(ints),function(i) list(z=c(),p1=c(.5)))
  
  counts <- ints[,7]
  d <- depth[,3]
  
  lambda0 <- c(1)
  lambda1 <- c(5)
  hsd0 <- hsd1 <- sd(counts)
  p1 <- c(.5)
  a0 <- c(1)
  b0 <- c(1)
  a1 <- c(1)
  b1 <- c(1)
  
  f <- function(p,h1,h0) p*dpois(counts,h1)/rowSums(cbind(p*dpois(counts,h1),(1-p)*dpois(counts,h0)))
  
  for( i in 1:N ) {
#    print(c(lambda1[i],lambda0[i]))
    vp <- f(sapply(l,function(v) v$p1[i]),lambda1[i],lambda0[i])
    print(vp)
    vz <- rbinom(S,1,vp)
    #print(sum(vz))
    l <- mapply(function(lx,z){
      p <- rbeta(1,1+z,1+(1-z))
      lx$z[i] <- z
      lx$p1[i+1] <- p
      lx
    },l,as.list(vz),SIMPLIFY=F)
    
    l0 <- 0
    r <- sum(counts[vz==0])
    #print(r)
   # print((a[i]*lambda0[i]+sum(vz==0)*r)/(a[i]+sum(vz)))
    while(l0 <= 0 ) l0 <- rgamma(1,a0[i]+r,b0[i]+sum(vz==0))
    a0[i+1] <- a0[i]+r
    b0[i+1] <- b0[i]+sum(vz==0)
#    print(l0)
    l1 <- l0
    r <- sum(counts[vz==1])
    #print(r)
    
    #print(c(a1[i]+r,b1[i]+sum(vz==1)))
    while(l1 <= l0 ) {
      l1 <- rgamma(1,a1[i]+r,b1[i]+sum(vz==1))

#      print(c(l0,l1))
    }
    a1[i+1] <- a1[i]+r
    b1[i+1] <- b1[i]+sum(vz==1)
    
    
    lambda0[i+1] <- l0
    lambda1[i+1] <- l1
    #print('---')
  }
  
  list(s=l,l0=lambda0,l1=lambda1,a0=a0,b0=b0,a1=a1,b1=b1)
}

extract.global.bayesian.prob <- function(model) {
  sapply(model$s,function(v) sum(v$z)/length(v$z))
}