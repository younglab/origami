library(GenomicRanges,quietly = !interactive())
library(utils,quietly = !interactive())
library(matrixStats,quietly = !interactive())

estimate.global.bayesian.mixture <- function(ints,depth,N=1100,burnin=100,pruning=NULL,with.distance.weight=F,multiply=T,show.progress=F,
                                             lambda0.init=1,lambda1.init=5) {
  S <- nrow(ints)
  
  d <- GRanges(seqnames=as.character(depth$V1),ranges=IRanges(depth$V2,depth$V3),strand='*')
  g1 <- GRanges(seqnames=as.character(ints$V1),ranges=IRanges(ints$V2,ints$V3),strand='*')
  g2 <- GRanges(seqnames=as.character(ints$V4),ranges=IRanges(ints$V5,ints$V6),strand='*')
  
  m1 <- match(g1,d)
  m2 <- match(g2,d)
  counts <- ints[,7]
  d <- depth[,4]
  
  if(!multiply) {
    sdepth <- rowSums(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  } else {
    sdepth <- rowProds(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  }
  
  print(c(msdepth,range(sdepth)))
  
  pp <- rep(.5,S)
  
  ret <- list(z=vector("list",N),
              p1=c(list(pp),vector("list",N)),
              mp=vector("list",N),
              lambda0=c(lambda0.init,rep(NA_real_,N)),
              lambda1=c(lambda1.init,rep(NA_real_,N)))

  
  totcounts <- sum(counts)

  
  if(show.progress) pb <- txtProgressBar()
  
  for( i in 1:N ) {
    lambda0 <- ret$lambda0[i]
    lambda1 <- ret$lambda1[i]
    g1 <- pp*dpois(counts,lambda1)
    g2 <- (1-pp) * dpois(counts,lambda0)
    
    vp <- g1/rowSums(cbind(g1,g2)) 

    if(any(is.na(vp))) {
      b<- is.na(vp)
      vp[b] <- ifelse(counts[b]>=lambda1,.999,.001)
    }
    
    vz <- rbinom(S,1,vp)

    pp <- rbeta(S,sdepth/msdepth+vz,1+(1-vz)) #msdepth/msdepth=1
    
    
    ret$z[[i]] <- vz
    ret$mp[[i]] <- vp
    ret$p1[[i+1]] <- pp
    
    r <- sum(counts[vz==0])
    n <- sum(vz==0)
    
    l0 <- rgamma(1,r,n)


    l1x <- rgamma(1,totcounts-r,S-n)#sum(vz==1))
    l1 <- max(l1,l1x)
    
    ret$lambda0[i+1] <- l0
    ret$lambda1[i+1] <- l1

    print(c(l0,l1))
    if(show.progress) setTxtProgressBar(pb,i/N)
  }
  
  if(show.progress) close(pb)
  if(!is.null(burnin) && is.numeric(burnin) && burnin > 0) {
    bn <- lapply(ret,function(l) l[1:burnin])
    idx <- -(1:burnin)
    ret <- lapply(ret,function(v) v[idx])
    ret$burnin <- bn
  }
  ret <- c(ret,list(sdepth=sdepth,msdepth=msdepth))
  ret
}

estimate.global.bayesian.mixture.candidate1 <- function(ints,depth,totint,N=1100,burnin=100,pruning=NULL,
                                                        with.distance.weight=F,multiply=T,show.progress=F,
                                                        lambda0.init=1,lambda1.init=5) {
  S <- nrow(ints)
  
  d <- GRanges(seqnames=as.character(depth$V1),ranges=IRanges(depth$V2,depth$V3),strand='*')
  g1 <- GRanges(seqnames=as.character(ints$V1),ranges=IRanges(ints$V2,ints$V3),strand='*')
  g2 <- GRanges(seqnames=as.character(ints$V4),ranges=IRanges(ints$V5,ints$V6),strand='*')
  
  m1 <- match(g1,d)
  m2 <- match(g2,d)
  counts <- ints[,7]
  d <- depth[,4]
  
  if(!multiply) {
    sdepth <- rowSums(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  } else {
    sdepth <- rowProds(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  }
  
  f1 <- paste(ints$V1,ints$V2,ints$V3,sep='_')
  f2 <- paste(ints$V4,ints$V5,ints$V6,sep='_')
  
  m1 <- match(f1,names(totint))
  m2 <- match(f2,names(totint))
  
  # frequency of interactions not including the present one
  nint <- rowMeans(cbind(totint[m1]-1,totint[m2]-1)) # since each interaction counts itself in the table in the parent script

  depthscore <- floor(msdepth/sdepth)
  alphaparam <- 1+depthscore+2*counts
  betaparam <- 1+nint
  
  #print(range(nint))
  #print(c(msdepth,range(sdepth)))
  
  pp <- rep(.5,S)
  
  ret <- list(z=vector("list",N),
              p1=c(list(pp),vector("list",N)),
              mp=vector("list",N),
              lambda0=c(lambda0.init,rep(NA_real_,N)),
              lambda1=c(lambda1.init,rep(NA_real_,N)))
  
  totcounts <- sum(counts)
  
  
  if(show.progress) pb <- txtProgressBar()
  #x <- sdepth/msdepth
  
  for( i in 1:N ) {
    lambda0 <- ret$lambda0[i]
    lambda1 <- ret$lambda1[i]
    g1 <- pp*dpois(counts,lambda1)
    g2 <- (1-pp) * dpois(counts,lambda0)
    
    vp <- g1/rowSums(cbind(g1,g2)) 
    
    # Sometimes the counts value in dpois is so extreme the loss of precision causes the value to be 0, need to correct for this
    # if in this case the counts value is closer to lambda1 than lambda0, give it a probability of .999, otherwise .001
    if(any(is.na(vp))) {
      b<- is.na(vp)
      vp[b] <- ifelse(counts[b]>=lambda1,.999,.001)
    }
    
    vz <- rbinom(S,1,vp)
    
    pp <- rbeta(S,alphaparam+vz,betaparam+(1-vz)) 
    
    
    ret$z[[i]] <- vz
    ret$mp[[i]] <- vp
    ret$p1[[i+1]] <- pp
    
    r0 <- sum(counts[vz==0])
    n <- sum(vz==0)
    
    l0 <- rgamma(1,r0,n)
    
    l1 <- l0
    r1 <- totcounts - r0 #sum(counts[vz==1])
    
    
    
    l1x <- rgamma(1,r1,S-n)#sum(vz==1))
    l1 <- max(l1,l1x)
    
    ret$lambda0[i+1] <- l0
    ret$lambda1[i+1] <- l1
    
    #print(c(r0,n,r1,S-n))
    #print(c(l0,l1))
    if(show.progress) setTxtProgressBar(pb,i/N)
  }
  
  if(show.progress) close(pb)
  #ret <- list(s=l,l0=lambda0,l1=lambda1)
  if(!is.null(burnin) && is.numeric(burnin) && burnin > 0) {
    bn <- lapply(ret,function(l) l[1:burnin])
    idx <- -(1:burnin)
    ret <- lapply(ret,function(v) v[idx])
    ret$burnin <- bn
  }
  ret <- c(ret,list(sdepth=sdepth,msdepth=msdepth))
  ret
}

estimate.global.bayesian.no.depth.mixture <- function(ints,depth,N=1100,burnin=100,pruning=NULL,with.distance.weight=F,multiply=T,show.progress=F) {
  S <- nrow(ints)
  
  d <- GRanges(seqnames=as.character(depth$V1),ranges=IRanges(depth$V2,depth$V3),strand='*')
  g1 <- GRanges(seqnames=as.character(ints$V1),ranges=IRanges(ints$V2,ints$V3),strand='*')
  g2 <- GRanges(seqnames=as.character(ints$V4),ranges=IRanges(ints$V5,ints$V6),strand='*')
  
  m1 <- match(g1,d)
  m2 <- match(g2,d)
  counts <- ints[,7]
  d <- depth[,4]
  
  
  if(!multiply) {
    sdepth <- rowSums(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  } else {
    sdepth <- rowProds(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  }
  
  
  print(c(msdepth,range(sdepth)))
  
  l <- lapply(1:S,function(i) list(z=rep(NA_integer_,N),p1=c(.5,rep(NA_real_,N)),mp=rep(NA_real_,N)))
  pp <- rep(.5,S)
  
  lambda0 <- c(1,rep(NA_real_,N))
  lambda1 <- c(5,rep(NA_real_,N))
  
  
  totcounts <- sum(counts)
  
  
  if(show.progress) pb <- txtProgressBar()
  
  
  for( i in 1:N ) {
    #vp <- f(pp,lambda1[i],lambda0[i])
    g1 <- pp*dpois(counts,lambda1[i])
    g2 <- (1-pp) * dpois(counts,lambda0[i])
    
    vp <- g1/rowSums(cbind(g1,g2)) 
    if(any(is.na(vp))) {
      b<- is.na(vp)
      vp[b] <- as.integer(counts[b]>=lambda1[i])
    }    
    vz <- rbinom(S,1,vp)
    
    pp <- rbeta(S,1+vz,1+(1-vz)) 
    
    l <- mapply(function(lx,z,mp,p){
      lx$z[i] <- z
      lx$mp[i] <- mp
      lx$p1[i+1] <- p
      lx
    },l,as.list(vz),as.list(vp),as.list(pp),SIMPLIFY=F)
    
    l0 <- 0
    r <- sum(counts[vz==0])
    n <- sum(vz==0)
    
    l0 <- rgamma(1,r,n)
    
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
  if(!is.null(burnin) && is.numeric(burnin) && burnin > 0) {
    orig <- ret
    idx <- -(1:burnin)
    ret <- lapply(ret,function(v) v[idx])
    ret$orig <- orig
  }
  ret <- c(ret,list(sdepth=sdepth,msdepth=msdepth))
  ret
}


extract.no.depth.bayesian.prob <- function(model) {
  sapply(model$s,function(v) sum(v$z)/length(v$z))
}

extract.global.bayesian.mixture.prob <- function(model) {
  m <- do.call(rbind,model$z)
  N <- nrow(m)
  colSums(m)/N
}
