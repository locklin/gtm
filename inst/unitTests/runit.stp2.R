### --- Test setup ---
library("gtm")

## tolerance
mytol = .Machine$double.eps^0.5

## original function; keeping around for cannon
gtm.stp2.old <-function (TT, noLatVarSmpl, noBasFn, s,kind="rect") {
    if (s <= 0) {
      stop("Argument s must have strict positive value")
    }
    gridXdim = sqrt(noLatVarSmpl)
    gridFIdim = sqrt(noBasFn)
    if ((gridXdim != floor(gridXdim)) || (gridFIdim != floor(gridFIdim))) {
      stop("Invalid number of basis functions or latent variable size.")
    }
    X = gtm.rctg(gridXdim, gridXdim)
    MU = gtm.rctg(gridFIdim, gridFIdim)
    MU = MU * (gridFIdim/(gridFIdim - 1))
    sigma = s * (MU[1, 1] - MU[2, 1])
    FI = gtm.gbf(MU, sigma, X)
    pciResult = gtm.pci.beta(TT, X, FI)
    list(X = X, MU = MU, FI = FI, W = pciResult$W, beta = pciResult$beta)
}

TT = matrix(3:61/20, ncol=1)
TT = cbind(TT, TT + 1.25 * sin(2*TT))
         
## setup and training
old = gtm.stp2.old(TT, 81, 25, 2)
new = gtm.stp2(TT, 81, 25, 2)

### --- Test functions ---


test.stp2.diff <- function() {  
  checkEqualsNumeric(old$X,new$X, tolerance= mytol)
  checkEqualsNumeric(old$beta,new$beta, tolerance= mytol)
  checkEqualsNumeric(old$MU,new$MU, tolerance= mytol)
  checkEqualsNumeric(old$FI,new$FI, tolerance= mytol)
  checkEqualsNumeric(old$W,new$W, tolerance= mytol)
}
