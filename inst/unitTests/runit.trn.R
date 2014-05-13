### --- Test setup ---
library("RUnit")
library("gtm")


## original function; keeping around for cannon
gtm.trn.old <-function(TT, FI, W, l, cycles, beta, m = 1, quiet = FALSE, minSing = 0.01) {
  loglik = matrix(0, cycles, 1)
  FI.T = t(FI)
  K = nrow(FI)
  Mplus1 = col(FI)
  N = nrow(TT)
  D = ncol(TT)
  ND = N * D
  A = matrix(0, Mplus1, Mplus1)
  cholDcmp = matrix(0, Mplus1, Mplus1)
  if (l > 0) {
    LAMBDA = l * matrix(1, Mplus1)
    LAMBDA[Mplus1] = 0
  }
  gtmGlobalDist = gtm.dist(TT, FI %*% W, m)
  if (m > 0) {
    gtmGlobalMinDist = gtmGlobalDist$minDist
    gtmGlobalMaxDist = gtmGlobalDist$maxDist
    gtmGlobalDist = gtmGlobalDist$DIST
  }
  for (cycle in 1:cycles) {
    llh = gtm.resp6(gtmGlobalDist, gtmGlobalMinDist, gtmGlobalMaxDist, 
      beta, D, m)
    gtmGlobalR = llh$R
    llh = llh$llh
    loglik[cycle] = llh
    if (!quiet) {
      matplot(1:cycle, loglik[1:cycle], xlim = c(1, cycles), 
              xlab = "Training cycle", ylab = "log-likelihood", 
              pch = 21)
      out1 = sprintf("Cycle: %d\tlogLH: %g\tBeta: %g\n", 
        cycle, llh, beta)
      cat(out1)
    }
    if (l > 0) {
      A = t(t(FI.T) * rowSums(gtmGlobalR)) %*% FI + LAMBDA/beta
    } else {
      A = t(t(FI.T) * rowSums(gtmGlobalR)) %*% FI
    }
    cholResult = chol(A, pivot = TRUE)
    if (attr(cholResult, "rank") != nrow(A)) {
      if (!quiet) {
        print("gtm.trn: Warning -- M-Step matrix singular, using pinv.\n")
      }
      svdResult = svd(A)
      sd = svdResult$d
      N = sum(ifelse(sd > minSing, 1, 0))
      warning(sprintf("Using %d out of %d eigenvalues", 
                      N, nrow(A)))
      if (N < 1) {
        stop("very singular matrix")
      }
      svdInverse = svdResult$v[, 1:N] %*% diag(1/sd[1:N], 
        nrow = N) %*% t(svdResult$u[, 1:N])
      W = svdInverse %*% (FI.T %*% (gtmGlobalR %*% TT))
    } else {
      oo = order(attr(cholResult, "pivot"))
      W = chol2inv(cholResult)[oo, oo] %*% (FI.T %*% (gtmGlobalR %*% TT))
    }
    gtmGlobalDist = gtm.dist(TT, FI %*% W, m)
    if (m > 0) {
      gtmGlobalMinDist = gtmGlobalDist$minDist
      gtmGlobalMaxDist = gtmGlobalDist$maxDist
      gtmGlobalDist = gtmGlobalDist$DIST
    }
    ##       gc() ## this takes a lot of time ... why is it there?
    beta = ND/sum(colSums(gtmGlobalDist * gtmGlobalR))
  }
  list(W = W, beta = beta, loglik = loglik)
}

TT = matrix(3:61/20, ncol=1)
TT = cbind(TT, TT + 1.25 * sin(2*TT))
         
## setup and training
res = gtm.stp1(TT, 20, 5, 2)
old = gtm.trn.old(TT, res$FI, res$W, 0.0, 20, res$beta,quiet=TRUE)
new = gtm.trn(TT, res$FI, res$W, 0.0, 20, res$beta,quiet=TRUE)

### --- Test functions ---


test.trn.diff <- function() {  
  checkEqualsNumeric(old$W,new$W,
                     tolerance= .Machine$double.eps^0.5)
    checkEqualsNumeric(old$beta,new$beta,
                     tolerance= .Machine$double.eps^0.5)
    checkEqualsNumeric(old$loglik,new$loglik,
                     tolerance= .Machine$double.eps^0.5)
}
