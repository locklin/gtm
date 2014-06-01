### --- Test setup ---
library("gtm")


## original function; keeping around for cannon
gtm.dist.old <-function (TT, Y, m = 0) {
  if (m < 0) {
    stop("Invalid value for mode ")
  }
  N = nrow(TT)
  tD = ncol(TT)
  K = nrow(Y)
  yD = ncol(Y)
  if (yD != tD) {
    stop("Mismatch in number of columns between T and Y.")
  }
  if (yD > 1) {
    DIST = Y %*% t(TT)
    ttt = matrix(colSums(t(TT)^2), nrow = 1)
    yy = matrix(colSums(t(Y)^2), nrow = 1)
    DIST = t(yy) %*% matrix(1, 1, N) + matrix(1, K, 1) %*% ttt - 2 * DIST
  } else {
    DIST = matrix(0, K, N)
    for (n in 1:N) {
      for (k in 1:K) {
        DIST[k, n] = sum((TT[n,] - Y[k, ])^2)
      }
    }
  }
  if (m > 0) {
    list(DIST = DIST, minDist = apply(DIST, 2, min), maxDist = apply(DIST, 2, max))
  } else {
    DIST
  }
}



A <- matrix(c(1,3,4,2,0,-1), ncol=2)
B <- matrix(c(1,0,0,1), ncol=2)

AA <- matrix(c(1,3,4,2,0,-1), ncol=1)
BB <- matrix(c(1,1,0,1), ncol=1)
 
### --- Test functions ---
 
test.proxy.dist.reg <- function() {
  checkEqualsNumeric(gtm.dist(A,B),gtm.dist.old(A,B),
                     tolerance=.Machine$double.eps^0.5)
}


test.proxy.dist.flat <- function(){
  checkEqualsNumeric(gtm.dist(AA,BB),gtm.dist.old(AA,BB),
                     tolerance=.Machine$double.eps^0.5)
}

test.proxy.dist.flat.str <- function(){
  gtm.dist.old(AA,BB,1) -> old
  gtm.dist(AA,BB,1) -> new
  checkEqualsNumeric(old$DIST,new$DIST,
                     tolerance=.Machine$double.eps^0.5)
  checkEqualsNumeric(old$minDist,new$minDist,
                     tolerance=.Machine$double.eps^0.5)
  checkEqualsNumeric(old$maxDist,new$maxDist,
                     tolerance=.Machine$double.eps^0.5)
}
  

test.proxy.dist.reg.str <- function() {
  gtm.dist.old(A,B,1) -> old
  gtm.dist(A,B,1) -> new
  checkEqualsNumeric(old$DIST,new$DIST,
                     tolerance=.Machine$double.eps^0.5)
  checkEqualsNumeric(old$minDist,new$minDist,
                     tolerance=.Machine$double.eps^0.5)
  checkEqualsNumeric(old$maxDist,new$maxDist,
                     tolerance=.Machine$double.eps^0.5)
}
