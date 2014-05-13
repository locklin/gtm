### --- Test setup ---
library("RUnit")
library("gtm")

A <- matrix(c(1,3,4,2,0,-1), ncol=2)
B <- matrix(c(1,0,0,1), ncol=2)

AA <- matrix(c(1,3,4,2,0,-1), ncol=1)
BB <- matrix(c(1,1,0,1), ncol=1)
 
### --- Test functions ---
 
test.proxy.dist.reg <- function() {
  checkTrue(gtm.dist(A,B)==gtm.dist.old(A,B))
}

test.proxy.dist.reg.str <- function() {
  checkTrue(gtm.dist(A,B,1)==gtm.dist.old(A,B,1))
}

test.proxy.dist.flat <- function(){
  checkTrue(gtm.dist(AA,BB)==gtm.dist.old(AA,BB))
}

test.proxy.dist.flat.str <- function(){
  old = gtm.dist.old(AA,BB,1)
  new = gtm.dist(AA,BB,1)
  checkTrue(old$DIST==new$DIST)
  checkTrue(old$minDist==new$minDist)
  checkTrue(old$maxDist==new$maxDist)
}
  

