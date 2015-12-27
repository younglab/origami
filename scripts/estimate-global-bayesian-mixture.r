estimate.global.bayesian.mixture <- function(ints,depth,N=1000) {
  S <- nrow(ints)
  
  l <- lapply(1:nrow(ints),function(i) list(z=c(),p1=c(.5)))
  
  counts <- ints[,7]
  d <- depth[,3]
  
  lambda0 <- c(1)
  lambda1 <- c(5)
  hsd0 <- hsd1 <- sd(counts)
  p1 <- c(.5)
  a <- c(1)
  b <- c(1)
  
  f <- function(p,h1,h0) p*dpois(counts,h1)/rowSums(cbind(p*dpois(counts,h1),(1-p)*dpois(counts,h0)))
  
  for( i in 1:N ) {
#    print(c(lambda1[i],lambda0[i]))
    vp <- f(sapply(l,function(v) v$p1[i]),lambda1[i],lambda0[i])
#    print(vp)
    vz <- rbinom(S,1,vp)
    print(sum(vz))
    l <- mapply(function(lx,z){
      p <- rbeta(1,1+z,1+(1-z))
      lx$z[i] <- z
      lx$p1[i+1] <- p
      lx
    },l,as.list(vz),SIMPLIFY=F)
    
    l0 <- 0
    r <- max(mean(counts[vz==0]),0,na.rm=T) ## default to 0
    #print(r)
   # print((a[i]*lambda0[i]+sum(vz==0)*r)/(a[i]+sum(vz)))
    while(l0 <= 0 ) l0 <- rnorm(1,(a[i]*lambda0[i]+sum(vz==0)*r)/(a[i]+sum(vz==0)),hsd0)
    a[i+1] <- a[i]+sum(vz==0)+1
    if(sum(vz==0) >= 2) hsd0 <- sd(c(counts[vz==0]))
#    print(l0)
    l1 <- l0
    r <- max(mean(counts[vz==1]),0,na.rm=T) ## default to 0
    
    while(l1 <= l0 ) {
      l1 <- rnorm(1,(b[i]*lambda1[i]+sum(vz==1)*r)/(b[i]+sum(vz==1)),hsd1)
      print(c(l0,l1))
    }
    b[i+1] <- b[i] + sum(vz==1)+1
    if( sum(vz==1) >= 2)hsd1 <- sd(c(counts[vz==1]))
   print(c(hsd0,hsd1))
    
    
    lambda0[i+1] <- l0
    lambda1[i+1] <- l1
    print('---')
  }
  
  list(s=l,l0=lambda0,l1=lambda1,a=a,b=b)
}

extract.global.bayesian.prob <- function(model) {
  sapply(model$s,function(v) sum(v$z)/length(v$z))
}