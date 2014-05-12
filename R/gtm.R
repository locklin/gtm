gtm.bi <-function (Y) {
    yInterDist = gtm.dist(Y, Y) + diag(10^10, nrow(Y))
    meanNN = mean(apply(yInterDist, 2, min))
    return(2/meanNN)
}


gtm.dist <-function (TT, Y, mode = 0) {
  if (mode < 0) {
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
  if (mode > 0) {
    return(list(DIST = DIST, minDist = apply(DIST, 2, min), 
                maxDist = apply(DIST, 2, max)))
  } else {
    return(DIST)
  }
}


## train a GTM
gtm.trn <-function (TT, FI, W, l, cycles, beta, m = 1, quiet = FALSE, minSing = 0.01) {
  mode = m
  llhLog = matrix(0, cycles, 1)
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
  gtmGlobalDist = gtm.dist(T, FI %*% W, mode)
  if (mode > 0) {
    gtmGlobalMinDist = gtmGlobalDist$minDist
    gtmGlobalMaxDist = gtmGlobalDist$maxDist
    gtmGlobalDist = gtmGlobalDist$DIST
  }
  for (cycle in 1:cycles) {
    llh = gtm.resp6(gtmGlobalDist, gtmGlobalMinDist, gtmGlobalMaxDist, 
      beta, D, mode)
    gtmGlobalR = llh$R
    llh = llh$llh
    llhLog[cycle] = llh
    if (!quiet) {
      matplot(1:cycle, llhLog[1:cycle], xlim = c(1, cycles), 
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
    gtmGlobalDist = gtm.dist(TT, FI %*% W, mode)
    if (mode > 0) {
      gtmGlobalMinDist = gtmGlobalDist$minDist
      gtmGlobalMaxDist = gtmGlobalDist$maxDist
      gtmGlobalDist = gtmGlobalDist$DIST
    }
    ##       gc() ## this takes a lot of time ... why is it there?
    beta = ND/sum(colSums(gtmGlobalDist * gtmGlobalR))
  }
  return(list(W = W, beta = beta, llhLog = llhLog))
}

gtm.gbf <-function (MU, sigma, X) {
    K = nrow(X)
    L = ncol(X)
    M = nrow(MU)
    L2 = ncol(MU)
    if (L != L2) 
        stop("Mismatch in dimensions of input argument matrices.")
    DIST = gtm.dist(MU, X)
    FI = exp((-1/(2 * sigma^2)) * DIST)
    return(cbind(FI, matrix(1, nrow = K, ncol = 1)))
}

gtm.hxg <-function (xDim, yDim) {
    if ((xDim < 2) || (yDim < 2) || (yDim != floor(yDim)) || 
        (xDim != floor(xDim))) {
      stop("Invalid grid dimensions")
    }
    r1 = 0:(xDim - 1)
    r2 = (yDim - 1):0
    fx = function(x, y) return(x)
    fry = function(x, y) return(y * sqrt(3)/2)
    X = matrix(outer(r1, r2, fx), ncol = yDim)
    Y = outer(r1, r2, fry)
    i = 2
    while (i <= yDim) {
        X[, i] = X[, i] + 0.5
        i = i + 2
    }
    grid = cbind(matrix(X, ncol = 1), matrix(Y, ncol = 1))
    maxVal = max(grid)
    grid = grid * (2/maxVal)
    maxXY = apply(grid, 2, max)
    grid[, 1] = grid[, 1] - maxXY[1]/2
    grid[, 2] = grid[, 2] - maxXY[2]/2
    return(grid)
}

gtm.lbf <-function (X) {
    return(cbind(X, matrix(1, nrow = nrow(X), ncol = 1)))
}

gtm.pci <-function (TT, X, FI) {
    K = nrow(X)
    L = ncol(X)
    Mplus1 = ncol(FI)
    if (K != nrow(FI)) {
      stop("wrong number of rows")
    }
    eV = eigen(cov(TT))
    A = eV$vectors[, 1:L] %*% diag(eV$values[1:L]^0.5, nrow = L)
    normX = (X - (matrix(1, K, L) %*% diag(colMeans(X), ncol(X)))) %*% 
        diag(1/apply(X, 2, sd), ncol(X))
    W = lsfit(FI, normX %*% t(A), intercept = FALSE)$coefficients
    W[Mplus1, ] = colMeans(TT)
    return(W)
}

gtm.pci.beta <-function (TT, X, FI) {
    K = nrow(X)
    L = ncol(X)
    Mplus1 = ncol(FI)
    if (K != nrow(FI)) {
      stop("wrong number of rows")
    }
    eV = eigen(cov(TT))
    A = eV$vectors[, 1:L] %*% diag(eV$values[1:L]^0.5, nrow = L)
    normX = (X - (matrix(1, K, L) %*% diag(colMeans(X), ncol(X)))) %*% 
        diag(1/apply(X, 2, sd), ncol(X))
    W = lsfit(FI, normX %*% t(A), intercept = FALSE)$coefficients
    W[Mplus1, ] = colMeans(TT)
    interDistBeta = gtm.bi(FI %*% W)
    if (L < length(TT[1, ])) {
      beta = min(interDistBeta, (1/eV$values[L + 1]))
    } else {
      beta = interDistBeta
    }
    return(list(W = W, beta = beta))
}

gtm.pmd <-function (TT, X, FI, W) {
    D = ncol(TT)
    DIST = gtm.dist(TT, FI %*% W)
    minDist = apply(DIST, 2, which.min)
    return(matrix(X[minDist, ], ncol = 1))
}

gtm.pmn <-function (TT, X, FI, W, b) {
    D = ncol(TT)
    DIST = gtm.dist(TT, FI %*% W)
    R = gtm.resp3(DIST, b, D)$R
    return(t(R) %*% X)
}

gtm.ppd1 <-function (tt, Y, beta, X) {
    tt = matrix(t, nrow = 1)
    D = ncol(tt)
    L = ncol(X)
    distResult = gtm.dist(tt, Y, 1)
    R = gtm.resp6(distResult$DIST, distResult$minDist, distResult$maxDist, 
        beta, D, 1)$R
    if (L == 1) {
      return(list(X = X, P = R))
    } else {
      stop("wrong argument")
    }
}

gtm.ppd2 <-function (tt, Y, beta, X, xDim, yDim) {
    tt = matrix(tt, nrow = 1)
    N = nrow(tt)
    D = ncol(tt)
    K = nrow(X)
    L = ncol(X)
    distResult = gtm.dist(matrix(tt, nrow = 1), Y, 1)
    R = gtm.resp6(distResult$DIST, distResult$minDist, distResult$maxDist, 
        beta, D, 1)$R
    if (L == 2) {
      return(gtm.r2m3(X[, 1], X[, 2], R, xDim, yDim))
    } else{
      stop("Wrong argument")
    }
}

gtm.pts <-function (M) {
    N = M - 1
    return(matrix(((-N/2):(N/2))/(N/2), ncol = 1))
}

gtm.r2m1 <-function (cX, meshRows, meshCols) {
    return(list(X = matrix(cX, nrow = meshRows)))
}

gtm.r2m2 <-function (cX, cY, meshRows, meshCols) {
    return(list(X = matrix(cX, nrow = meshRows), Y = matrix(cY, 
        nrow = meshRows)))
}

gtm.r2m3 <-function (cX, cY, cZ, meshRows, meshCols) {
    return(list(X = matrix(cX, nrow = meshRows), Y = matrix(cY, 
        nrow = meshRows), Z = matrix(cZ, nrow = meshRows)))
}

gtm.rctg <-function (xDim, yDim) {
    if ((xDim < 2) || (yDim < 2) || (yDim != floor(yDim)) || 
        (xDim != floor(xDim))) {
        stop("Invalid grid dimensions")
      }
    r1 = 0:(xDim - 1)
    r2 = (yDim - 1):0
    fx = function(x, y) return(x)
    fy = function(x, y) return(y)
    X = outer(r1, r2, fx)
    Y = outer(r1, r2, fy)
    grid = cbind(matrix(X, ncol = 1), matrix(Y, ncol = 1))
    maxVal = max(grid)
    grid = grid * (2/maxVal)
    maxXY = apply(grid, 2, max)
    grid[, 1] = grid[, 1] - maxXY[1]/2
    grid[, 2] = grid[, 2] - maxXY[2]/2
    return(grid)
}

gtm.resp3 <- function (DIST, beta, D) {
    return(gtm.resp6(DIST, beta, D, beta, D, 0))
}

gtm.resp6 <- function (DIST, minDist, maxDist, beta, D, mode) {
    if (mode < 0 || mode > 2 || mode != floor(mode)) {
      stop("Unknown mode of calculation")
    }
    if (is.matrix(beta) || is.matrix(D)) {
      stop("beta and D should be scalars - mismatch of arguments?")
    }
    if (D != floor(D)) {
      stop("Invalid value for D ")
    }
    K = nrow(DIST)
    N = ncol(DIST)
    if (mode > 0) {
        distCorr = (maxDist + minDist)/2
        distCorr = pmin(distCorr, (minDist + 700 * (2/beta)))
        for (n in 1:N) {
          DIST[, n] = DIST[, n] - distCorr[n]
        }
    }
    R = exp((-beta/2) * DIST)
    if (mode < 2) {
      rSum = colSums(R)
    } else {
      rSum = colSums(gtm.sort(R))
    }
    for (n in 1:N) {
      R[, n] = R[, n]/rSum[n]
    }
    if (mode < 1) {
      llh = sum(log(rSum)) + N * ((D/2) * log(beta/(2 * pi)) - log(K))
    } else {
      llh = sum(log(rSum) + distCorr * (-beta/2)) + N * ((D/2) * 
        log(beta/(2 * pi)) - log(K))
    }
    return(list(llh = llh, R = R))
}

gtm.ri <-function (TT, FI) {
    N = nrow(TT)
    D = ncol(TT)
    K = nrow(FI)
    Mplus1 = ncol(FI)
    varT = matrix(apply(TT, 2, sd)^2, 1, D)
    mnVarFI = mean(apply(FI[, 1:(Mplus1 - 1)], 2, sd)^2)
    stdW = varT/(mnVarFI * (Mplus1 - 2))
    W = rbind(matrix(rnorm((Mplus1 - 1) * D), Mplus1 - 1, D) %*% 
        diag(c(sqrt(stdW)), D), matrix(0, 1, D))
    W[Mplus1, ] = colMeans(TT) - colMeans(FI %*% W)
    return(W)
}

gtm.sort <-function (R) {
    idx = matrix(0, 1, ncol(R))
    for (n in 1:ncol(R)) {
      idx[n] = sum(R[, n]^2)
    }
    return(R[, order(idx)])
}

gtm.stp1 <-function (TT, noLatVarSmpl, noBasFn, s) {
    if (floor(noLatVarSmpl) != noLatVarSmpl || floor(noBasFn) != 
        noBasFn || noLatVarSmpl < 0 || noBasFn < 0) {
      stop("Incorrect arguments")
    }
    if (s <= 0) {
      stop("Argument s must have strict positive value")
    }
    X = gtm.pts(noLatVarSmpl)
    MU = gtm.pts(noBasFn)
    MU = MU * (noBasFn/(noBasFn - 1))
    sigma = s * (MU[2] - MU[1])
    FI = gtm.gbf(MU, sigma, X)
    pciResult = gtm.pci.beta(TT, X, FI)
    return(list(X = X, MU = MU, FI = FI, W = pciResult$W, beta = pciResult$beta))
}

gtm.stp2 <-function (TT, noLatVarSmpl, noBasFn, s) {
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
    return(list(X = X, MU = MU, FI = FI, W = pciResult$W, beta = pciResult$beta))
}


## demo
gtm.demo <-function (randomSetup = FALSE, points = 5, samples = 20) {
    cat("\nYou've started the GTM demo, please wait while ")
    cat("data is being generated.\n\n")
    TT = matrix(3:61/20, ncol = 1)
    TT = cbind(TT, TT + 1.25 * sin(2 * TT))
    plot(TT[, 1], TT[, 2], col = "red", pch = 21, xlim = c(0, 3.5), 
        ylim = c(0, 3.5))
    cat("TThe figure shows data generated by feeding a 1D uniform distribution\n")
    cat("(on the X-axis) through a non-linear function (y = x + 1.25*sin(2*x))\n")
    cat("\nPress Enter to continue ...\n\n")
    readString = readline()
    cat("Please wait while the GTTM model is set up.\n\n")
    if (randomSetup) {
        X = gtm.pts(samples)
        MU = gtm.pts(points)
        sigma = 1
        FI = gtm.gbf(MU, sigma, X)
        W = gtm.ri(TT, FI)
        b = gtm.bi(Y)
    } else {
        stpResult = gtm.stp1(TT, samples, points, 2)
        b = stpResult$beta
        FI = stpResult$FI
        W = stpResult$W
        X = stpResult$X
    }

    Y = FI %*% W
    matplot(Y[, 1], Y[, 2], pch = 23, col = "green", type = "l", 
        add = TRUE)
    plot(TT[, 1], TT[, 2], col = "red", pch = 21, xlim = c(0, 3.5), 
        ylim = c(0, 3.5))
    matplot(Y[, 1], Y[, 2], pch = 23, col = "green", type = "p", 
        add = TRUE)
    matplot(Y[, 1], Y[, 2], col = "green", type = "l", add = TTRUE)
    symbols(x = Y[, 1], y = Y[, 2], circles = matrix(sqrt(1/b), 
        ncol = 1, nrow = nrow(Y)), fg = "red", add = TRUE, inches = FALSE)
    title("Initial configuration")
    cat("The figure shows the starting point for the GTM, before the training.\n")
    cat("A discrete latent variable distribution of 20 points in 1 dimension \n")
    cat("is mapped to the 1st principal component of the target data.\n")
    cat("Each of the 20 points defines the centre of a Gaussian in a Gaussian \n")
    cat("mixture, marked by the green '+'-signs. The mixture components have \n")
    cat("all equal variance, illustrated by the filled circle around each \n")
    cat("'+'-sign, the raddii corresponding to 2 standard deviations.\n")
    cat("The '+'-signs are connected with a line according to their \n")
    cat("corresponding ordering in latent space.\n\n")
    cat("Press any key to begin training ...\n\n")
    cat("\nPress Enter to continue ...\n\n")
    readString = readline()
    for (j in 1:15) {
        trnResult = gtm.trn(TT, FI, W, 0, 1, b, quiet = TTRUE)
        W = trnResult$W
        b = trnResult$beta
        Y = FI %*% W
        plot(TT[, 1], TT[, 2], col = "red", xlim = c(0, 3.5), ylim = c(0, 
            3.5))
        matplot(Y[, 1], Y[, 2], col = "green", type = "p", add = TTRUE, 
            pch = 23)
        matplot(Y[, 1], Y[, 2], col = "green", type = "l", add = TTRUE)
        symbols(x = Y[, 1], y = Y[, 2], circles = matrix(sqrt(1/b), 
            ncol = 1, nrow = nrow(Y)), fg = "red", add = TTRUE, 
            inches = FALSE)
        title(sprintf("After %d iterations of training.", j))
        if (j == 4) {
            cat("TThe GTM initiallaly adapts relatively quickly - already after \n")
            cat("4 iterations of training, a rough fit is attained.\n\n")
            cat("\nPress Enter to continue ...\n\n")
            readString = readline()
        } else if (j == 8) {
            cat("After another 4 iterations of training:  from now on further \n")
            cat("training only makes small changes to the mapping, which combined with \n")
            cat("decrements of the Gaussian mixture variance, optimize the fit in \n")
            cat("terms of likelihood.\n\n")
            cat("Press any key to continue training ...\n\n")
            cat("\nPress Enter to continue ...\n\n")
            readString = readline()
        } else {
          Sys.sleep(1)
        }
    }
    cat("After 15 iterations of training the GTM can be regarded as converged. \n")
    cat("Is has been adapted to fit the target data distribution as well \n")
    cat("as possible, given prior smoothness constraints on the mapping. It \n")
    cat("captures the fact that the probabilty density is higher at the two \n")
    cat("bends of the curve, and lower towards its end points.\n\n")
    cat("Thanks for your interest in running the GTM demo.\n\n")
    return(list(X = X, TT = TT, W = W, Y = Y, beta = b, llh = trnResult$llh, FI = FI))
}
