library(GenomicRanges,quietly = !interactive())
library(utils,quietly = !interactive())
library(matrixStats,quietly = !interactive())
library(MASS)
library(mvtnorm)

calc.zscore <- function(v) (v-mean(v))/sd(v)

estimate.global.bayesian.mixture <- function(ints,depth,inttable,N=1100,burnin=100,pruning=NULL,
                                                        with.distance.weight=F,multiply=T,show.progress=F,
                                                        lambda0.init=1,lambda1.init=5,suppress.counts.higher.than=30,
                                                        mini.model=T,usedf=0) {
  S <- nrow(ints)
  
  d <- GRanges(seqnames=as.character(depth$V1),ranges=IRanges(depth$V2,depth$V3),strand='*')
  g1 <- GRanges(seqnames=as.character(ints$V1),ranges=IRanges(ints$V2,ints$V3),strand='*')
  g2 <- GRanges(seqnames=as.character(ints$V4),ranges=IRanges(ints$V5,ints$V6),strand='*')
  
  m1 <- match(g1,d)
  m2 <- match(g2,d)
  counts <- ints[,7]
  d <- depth[,4]
  intdist <- distance(g1,g2)
  interchromosomal <- is.na(intdist)
  minintdist <- min(intdist[!interchromosomal])
  
  pp1avg <- rep(0,S)
  pp1avgN <- 0
  
  proposal.var <- 10
  dbeta1.acs <- dbeta0.acs <- 0
  
  initial.beta <- list(matrix(c(0,0),ncol=1))

  
  if(is.null(pruning) || pruning < 1) {
    pruning <- N+1 # will never prune
  } 
  
  if(!multiply) {
    sdepth <- rowSums(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  } else {
    sdepth <- rowProds(cbind(d[m1],d[m2]))
    msdepth <- median(sdepth)
  }
  
  f1 <- paste(ints$V1,ints$V2,ints$V3,sep='_')
  f2 <- paste(ints$V4,ints$V5,ints$V6,sep='_')
  
  m1 <- match(f1,rownames(inttable))
  m2 <- match(f2,rownames(inttable))
  
  # frequency of interactions not including the present one
  #nint <- rowMeans(cbind(totint[m1]-1,totint[m2]-1)) # since each interaction counts itself in the table in the parent script
  
  #pint <- apply(cbind(inttable[m1],inttable[m2],counts),1,function(r) { idx <- as.integer(colnames(inttable))<r[3]; sum(r[1:2])}
  countm <- t(mapply(function(x1,x2,n) {
    idx <- as.integer(colnames(inttable))<n
    
    z <- inttable[c(x1,x2),]

    a <- if(!any(idx)) 0 else sum(z[,idx])
    b <- sum(z)-a
    c(a,b)
  },as.list(m1),as.list(m2),as.list(counts),SIMPLIFY=T))
  
  
  depthscore <- floor(sdepth/msdepth)
  alphaparam <- 1+countm[,1]
  betaparam <- 1+depthscore+countm[,2]
  

  
  pp <- rep(.5,S)
  
  if (!mini.model) {
    ret <- list(
      z = vector("list",N),
      p1 = c(list(pp),vector("list",N)),
      mp = vector("list",N),
      lambda0 = c(lambda0.init,rep(NA_real_,N)),
      lambda1 = c(lambda1.init,rep(NA_real_,N)),
      dbeta0 = c(initial.beta,vector("list",N)),
      dbeta1 = c(initial.beta,vector("list",N)),
      lambdad1 = vector("list",N),
      lambdad0 = vector("list",N)
    )
  } else {
    ret <- list(
      lambda0 = c(lambda0.init,rep(NA_real_,N)),
      lambda1 = c(lambda1.init,rep(NA_real_,N)),
      dbeta0 = c(initial.beta,vector("list",N)),
      dbeta1 = c(initial.beta,vector("list",N))
    )
    zm <- rep(0,length(counts))
  }
  
  totcounts <- sum(counts)
  
  lambdad1 <- lambdad0 <- rep(0,length(counts))
  
  
  if(show.progress) pb <- txtProgressBar()
  
  suppress <- counts > suppress.counts.higher.than
  
  
  for( i in 1:N ) {
    lambda0 <- ret$lambda0[i]
    lambda1 <- ret$lambda1[i]
    g1 <- pp*dpois(counts,lambda1 + lambdad1)
    g2 <- (1-pp) * dpois(counts,lambda0 + lambdad0)
    
    vp <- g1/rowSums(cbind(g1,g2)) 
    
    # Sometimes the counts value in dpois is so extreme the loss of precision causes the value to be 0, need to correct for this
    # if in this case the counts value is closer to lambda1 than lambda0, give it a probability of .999, otherwise .001
    if(any(is.na(vp))) {
      b<- is.na(vp)
      vp[b] <- ifelse(counts[b]>=lambda1,.999,.001)
    }
    
    vz <- rbinom(S,1,vp)
    
    pp <- rbeta(S,alphaparam+vz,betaparam+(1-vz)) 
    
    
    if(!mini.model) {
      ret$z[[i]] <- vz
      ret$mp[[i]] <- vp
      ret$p1[[i+1]] <- pp
    } else if(mini.model && i > burnin && ((i-burnin) %% pruning) != 0 ) {
      zm <- zm + vz
      pp1avg <- (pp1avg*pp1avgN+pp)/(pp1avgN+1)
      pp1avgN <- pp1avgN+1
    }
    
    b <- vz == 0 & !suppress
    r0 <- sum(counts[b])
    n <- sum(b)

    l0 <- rgamma(1,r0,n)
    
    l1 <- l0
    
    b <- vz == 1 & !suppress
    r1 <- sum(counts[b])
    n <- sum(b)
    
    
    
    l1x <- rgamma(1,r1,n)
    l1 <- max(l1,l1x)
    
    ret$lambda0[i+1] <- l0
    ret$lambda1[i+1] <- l1
    
    if(with.distance.weight) {
      dbeta1 <- ret$dbeta1[[i]]
      dbeta0 <- ret$dbeta0[[i]]
      
      d <- log10(intdist[vz==1 & !interchromosomal]+1)
      idx <- sample.int(length(d),min(1000,length(d)))
      x <- cbind(rep(1,length(d)),d)
      
      pc <- floor(pmax(counts[vz==1&!is.na(intdist)],0))
      
      #print(dbeta1)
      #print(sd(log(pc[idx]+.5)) *solve(t(x)%*%x))
      dbeta1.p <- t(rmvnorm(1,dbeta1,sd(log(pc[idx]+.5)) *solve(t(x)%*%x)))#*((x*t(x))^-1)[idx])
      #print(head(pc))
      #print(head(x))
      #print(head(t((x*t(x))^-1)))
      #print(head(t(var(log(pc+.5))*(x*t(x))^-1),100))
      #cat(paste("dbeta1s",dbeta1,dbeta1.p),sep='\n')
      #print(dbeta1)
      #print(dbeta1.p)
      
      lhr <- sum(dpois(pc,exp(x%*%dbeta1.p),log = T)) - sum(dpois(pc,exp(x%*%dbeta1),log=T)) + 
        sum(dnorm(dbeta1.p,rep(0,length(dbeta1)),rep(10,length(dbeta1)),log=T)) - sum(dnorm(dbeta1,rep(0,length(dbeta1)),rep(10,length(dbeta1)),log=T))
      #print(head(cbind(pc,x,dpois(pc,exp(dbeta1.p*x),log = T),1000)))
      #print(head(cbind(pc,x,dpois(pc,exp(dbeta1*x),log=T),1000)))
      #print(lhr)
      if(log(runif(1))<lhr) { dbeta1 <- dbeta1.p; dbeta1.acs <- dbeta1.acs+1 }
      cat(paste("dbeta1s",dbeta1,dbeta1.p,lhr),sep='\n')
      
      #s1 <- if( usedf > 0 ) smooth.spline(x,pmax(counts[vz==1& !is.na(intdist)]-l1,0),df=usedf) else smooth.spline(x,pmax(counts[vz==1& !is.na(intdist)]-l1,0))
      d <- log10(intdist[vz==0 & !interchromosomal]+1)
      idx <- sample.int(length(d),min(1000,length(d)))
      x <- cbind(rep(1,length(d)),d)
      
      
      
      pc <- floor(pmax(counts[vz==0&!is.na(intdist)],0))
      
      dbeta0.p <- t(rmvnorm(1,dbeta0,sd(log(pc[idx]+.5)) *solve(t(x)%*%x)))
      #print(dbeta0)
      #print(dbeta0.p)
      

      
      lhr <- sum(dpois(pc,exp(x%*%dbeta0.p),log = T)) - sum(dpois(pc,exp(x%*%dbeta0),log=T)) + 
        sum(dnorm(dbeta0.p,rep(0,length(dbeta0)),rep(10,length(dbeta0)),log=T)) - sum(dnorm(dbeta0,rep(0,length(dbeta0)),rep(10,length(dbeta0)),log=T))
      
      if(log(runif(1))<lhr) { dbeta0 <- dbeta0.p; dbeta0.acs <- dbeta0.acs+1 }
      cat(paste("dbeta0s",dbeta0,dbeta0.p,lhr),sep='\n')
      
      ret$dbeta1[[i+1]] <- dbeta1
      ret$dbeta0[[i+1]] <- dbeta0
      
      #s0 <- if( usedf > 0 ) smooth.spline(x,pmax(counts[vz==0& !is.na(intdist)]-l0,0),df=usedf) else smooth.spline(x,pmax(counts[vz==0& !is.na(intdist)]-l0,0))
      
      d <- log10(intdist+1)
      if(any(interchromosomal)) d[interchromosomal] <- log10(minintdist+1) ## set eact interchromsomal interaction to shortest distance (which should have the highest mean read count)
      x <- cbind(rep(1,length(intdist)),d)
      
      #lambdad1 <- pmax(predict(s1,x)$y,0) ### floor the value at 0
      #lambdad0 <- pmax(predict(s0,x)$y,0) 
      
      lambdad1 <- rpois(length(x),exp(x %*% dbeta1))
      lambdad0 <- rpois(length(x),exp(x %*% dbeta0))
      

      
      if( !mini.model) {
        ret$lambdad1[[i]] <- lambdad1
        ret$lambdad0[[i]] <- lambdad0
      }
    }
    if(show.progress) setTxtProgressBar(pb,i/N)
  }
  
  if(show.progress) close(pb)
  #ret <- list(s=l,l0=lambda0,l1=lambda1)
  if(!is.null(burnin) && is.numeric(burnin) && burnin > 0 && burnin < N) {
      bn <- lapply(ret,function(l) l[1:burnin])
      idx <- -(1:burnin)
      ret <- lapply(ret,function(v) v[idx])
      ret$burnin <- bn
  }
  if(!is.null(pruning) && pruning > 1) {
    l <- length(ret$lambda0)
    
    idx <- seq(1,l,by=pruning)
    w <- which(names(ret)=="burnin")
    if(length(w)>0) {
      bn <- ret$burnin
      ret <- lapply(ret[-w],function(l) l[-idx])
      ret$burnin <- bn
    } else {
      ret <- lapply(ret,function(l) l[-idx])
    }
  } 
  ret <- c(ret,list(sdepth=sdepth,msdepth=msdepth,intdist=intdist,dbeta1.acs=dbeta1.acs,dbeta0.acs=dbeta0.acs))
  if(mini.model) ret <- c(ret,list(zm=zm,lambdad1=lambdad1,lambdad0=lambdad0,pp1avg=pp1avg))
  
  ret
}

extract.global.bayesian.mixture.prob <- function(model) {
  if( "zm" %in% names(model)) {
    return(model$zm / length(model$lambda0))
  } else {
    m <- do.call(rbind,model$z)
    N <- nrow(m)
    return(colSums(m)/N)
  }
}
