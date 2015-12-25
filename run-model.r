convert.to.factor <- function(pos,idx,l=NULL) {
  d <- pos[,idx]
  s <- paste(d[,1],d[,2],d[,3],sep='_')
  
  r <- if(is.null(l)) factor(s) else factor(s,levels=levels(l)) 
  r
}


depth <- read.table("peak-counts.txt",sep='\t')
intcounts <- read.table("int-counts.txt",sep='\t')

f <- convert.to.factor(depth,1:3)

i1 <- convert.to.factor(intcounts,1:3,f)
i2 <- convert.to.factor(intcounts,4:6,f)

p <- intcounts[i1 != i2,]

params <- lapply(1:nrow(p),function(i) list(I=c(0),w=c(.5),lambda1=2))

N <- 1000

f <- function(x,i,w) {
  if(i==1 ) w * dpois(x,lambda=5) else (1-w) * dpois(x,lambda=1)
}

runupdate <- function(v,l) {
  n <- length(l[[1]])
  
  p <- l$w[[n]]*dpois(v,lambda=l$lambda1[[n]])/sum(l$w[[n]]*dpois(v,lambda=l$lambda1[[n]])+(1-l$w[[1]])*dpois(v,lambda=1))
  #print(p)
  i <- rbinom(1,1,p)
  w <- rbeta(1,1+i,1+(1-i))
  lambda1 <- 1
  while(lambda1 <= 1) lambda1 <- rnorm(1,weighted.mean(c(1,l$lambda1[[n]]),c(1-w,w)),1)
  
  l$I[[n+1]] <- i
  l$w[[n+1]] <- w
  l$lambda1[[n+1]] <- lambda1
  
  l
}

for( i in 1:N) {
  print(i)
  for( j in 1:50) {
    params[[j]] <- runupdate(p[j,7],params[[j]])
  }
}
