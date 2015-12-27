estimate.per.sample.bayesian.probability <- function(ints,depth,N=1000) {
  params <-
    lapply(1:nrow(p),function(i)
      list(
        I = c(0),w = c(.5),lambda1 = 2
      ))
  
  runupdate <- function(v,l) {
    n <- length(l[[1]])
    
    p <-
      l$w[[n]] * dpois(v,lambda = l$lambda1[[n]]) / sum(l$w[[n]] * dpois(v,lambda =
                                                                           l$lambda1[[n]]) + (1 - l$w[[1]]) * dpois(v,lambda = 1))
    #print(p)
    i <- rbinom(1,1,p)
    w <- rbeta(1,1 + i,1 + (1 - i))
    lambda1 <- 1
    while (lambda1 <= 1)
      lambda1 <- rnorm(1,weighted.mean(c(1,l$lambda1[[n]]),c(1 - w,w)),1)
    
    l$I[[n + 1]] <- i
    l$w[[n + 1]] <- w
    l$lambda1[[n + 1]] <- lambda1
    
    l
  }
  
  for (i in 1:N) {
    print(i)
    for (j in 1:nrow(ints)) {
      params[[j]] <- runupdate(ints[j,7],params[[j]])
    }
  }
  
  sapply(params,function(l) sum(l$I)/length(l$I))
}
