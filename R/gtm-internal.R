# Internal functions for the gtm package
# All functions here are unexported; used by GtmFit, predict.gtm, plot.gtm

#' Squared Euclidean Distance Matrix
#'
#' Computes squared Euclidean distance matrix between rows of two matrices
#' using the proxy package.
#'
#' @param x Numeric matrix (N1 x D).
#' @param y Numeric matrix (N2 x D).
#' @param m Integer mode: 0 = return matrix only; > 0 = return list with
#'   DIST, minDist, maxDist.
#' @param kind Distance method string for proxy::dist. Default "Euclidean".
#' @return If m == 0, a K x N distance matrix. If m > 0, a list with
#'   elements DIST, minDist, maxDist.
#' @importFrom proxy dist
#' @keywords internal
.gtm_dist <- function(x, y, m = 0, kind = "Euclidean") {
  out <- t(proxy::dist(x, y, method = kind)^2)
  if (m > 0) {
    list(DIST = out, minDist = apply(out, 2, min),
         maxDist = apply(out, 2, max))
  } else {
    out
  }
}

#' Generate 1D Latent Grid Points
#'
#' Creates M evenly spaced points in \[-1, 1\].
#'
#' @param M Integer, number of grid points.
#' @return Numeric matrix (M x 1).
#' @keywords internal
.gtm_pts <- function(M) {
  N <- M - 1
  matrix(((-N / 2):(N / 2)) / (N / 2), ncol = 1)
}

#' Generate Rectangular Grid
#'
#' Creates a 2D rectangular grid centered at the origin, scaled to \[-1, 1\].
#'
#' @param xdim Integer, number of columns.
#' @param ydim Integer, number of rows.
#' @return Numeric matrix (xdim*ydim x 2).
#' @keywords internal
.gtm_rctg <- function(xdim, ydim) {
  if (xdim < 2 || ydim < 2 || ydim != floor(ydim) || xdim != floor(xdim)) {
    stop("Invalid grid dimensions")
  }
  r1 <- 0:(xdim - 1)
  r2 <- (ydim - 1):0
  fx <- function(x, y) x
  fy <- function(x, y) y
  X <- outer(r1, r2, fx)
  Y <- outer(r1, r2, fy)
  grid <- cbind(as.vector(X), as.vector(Y))
  max_val <- max(grid)
  grid <- grid * (2 / max_val)
  max_xy <- apply(grid, 2, max)
  grid[, 1] <- grid[, 1] - max_xy[1] / 2
  grid[, 2] <- grid[, 2] - max_xy[2] / 2
  grid
}

#' Generate Hexagonal Grid
#'
#' Creates a 2D hexagonal grid centered at the origin, scaled to \[-1, 1\].
#'
#' @param xdim Integer, number of columns.
#' @param ydim Integer, number of rows.
#' @return Numeric matrix (xdim*ydim x 2).
#' @keywords internal
.gtm_hex <- function(xdim, ydim) {
  if (xdim < 2 || ydim < 2 || ydim != floor(ydim) || xdim != floor(xdim)) {
    stop("Invalid grid dimensions")
  }
  r1 <- 0:(xdim - 1)
  r2 <- (ydim - 1):0
  fx <- function(x, y) x
  fry <- function(x, y) y * sqrt(3) / 2
  xx <- matrix(outer(r1, r2, fx), ncol = ydim)
  yy <- outer(r1, r2, fry)
  i <- 2
  while (i <= ydim) {
    xx[, i] <- xx[, i] + 0.5
    i <- i + 2
  }
  grid <- cbind(as.vector(xx), as.vector(yy))
  max_val <- max(grid)
  grid <- grid * (2 / max_val)
  max_xy <- apply(grid, 2, max)
  grid[, 1] <- grid[, 1] - max_xy[1] / 2
  grid[, 2] <- grid[, 2] - max_xy[2] / 2
  grid
}

#' Gaussian Basis Functions
#'
#' Evaluates Gaussian RBFs centered at MU on grid points X, with a bias column.
#'
#' @param MU Numeric matrix (M x L), basis function centers.
#' @param sigma Numeric > 0, basis width.
#' @param X Numeric matrix (K x L), grid points.
#' @return Numeric matrix (K x (M+1)), with bias column appended.
#' @keywords internal
.gtm_gbf <- function(MU, sigma, X) {
  K <- nrow(X)
  L <- ncol(X)
  M <- nrow(MU)
  L2 <- ncol(MU)
  if (L != L2) {
    stop("Mismatch in dimensions of input argument matrices.")
  }
  DIST <- .gtm_dist(MU, X)
  FI <- exp((-1 / (2 * sigma^2)) * DIST)
  cbind(FI, matrix(1, nrow = K, ncol = 1))
}

#' Linear Basis Functions
#'
#' Returns X with a bias column appended.
#'
#' @param X Numeric matrix (K x L).
#' @return Numeric matrix (K x (L+1)).
#' @keywords internal
.gtm_lbf <- function(X) {
  cbind(X, matrix(1, nrow = nrow(X), ncol = 1))
}

#' PCA Weight Initialization
#'
#' Initializes weight matrix W so that the mapped points span the first L
#' principal components of the data.
#'
#' @param TT Numeric matrix (N x D), training data.
#' @param X Numeric matrix (K x L), latent grid.
#' @param FI Numeric matrix (K x (M+1)), basis function activations.
#' @return Numeric matrix ((M+1) x D), weight matrix.
#' @importFrom stats cov lsfit
#' @keywords internal
.gtm_pci <- function(TT, X, FI) {
  K <- nrow(X)
  L <- ncol(X)
  mplus <- ncol(FI)
  if (K != nrow(FI)) {
    stop("wrong number of rows")
  }
  eV <- eigen(cov(TT))
  A <- eV$vectors[, 1:L, drop = FALSE] %*%
    diag(eV$values[1:L]^0.5, nrow = L)
  norm_X <- (X - matrix(1, K, L) %*% diag(colMeans(X), ncol(X))) %*%
    diag(1 / apply(X, 2, sd), ncol(X))
  W <- lsfit(FI, norm_X %*% t(A), intercept = FALSE)$coefficients
  W[mplus, ] <- colMeans(TT)
  W
}

#' PCA Weight + Beta Initialization
#'
#' Like `.gtm_pci` but also initializes the inverse variance beta.
#'
#' @param TT Numeric matrix (N x D), training data.
#' @param X Numeric matrix (K x L), latent grid.
#' @param FI Numeric matrix (K x (M+1)), basis function activations.
#' @return List with W (weight matrix) and beta (scalar).
#' @importFrom stats cov lsfit
#' @keywords internal
.gtm_pci_beta <- function(TT, X, FI) {
  K <- nrow(X)
  L <- ncol(X)
  mplus <- ncol(FI)
  if (K != nrow(FI)) {
    stop("wrong number of rows")
  }
  eV <- eigen(cov(TT))
  A <- eV$vectors[, 1:L, drop = FALSE] %*%
    diag(eV$values[1:L]^0.5, nrow = L)
  norm_X <- (X - matrix(1, K, L) %*% diag(colMeans(X), ncol(X))) %*%
    diag(1 / apply(X, 2, sd), ncol(X))
  W <- lsfit(FI, norm_X %*% t(A), intercept = FALSE)$coefficients
  W[mplus, ] <- colMeans(TT)
  inter_dist_beta <- .gtm_bi(FI %*% W)
  if (L < ncol(TT)) {
    beta <- min(inter_dist_beta, 1 / eV$values[L + 1])
  } else {
    beta <- inter_dist_beta
  }
  list(W = W, beta = beta)
}

#' Random Weight Initialization
#'
#' Initializes weight matrix W randomly, scaled to match data variance.
#'
#' @param TT Numeric matrix (N x D), training data.
#' @param FI Numeric matrix (K x (M+1)), basis function activations.
#' @return Numeric matrix ((M+1) x D).
#' @importFrom stats rnorm sd
#' @keywords internal
.gtm_ri <- function(TT, FI) {
  N <- nrow(TT)
  D <- ncol(TT)
  K <- nrow(FI)
  mplus <- ncol(FI)
  var_T <- matrix(apply(TT, 2, sd)^2, 1, D)
  mn_var_FI <- mean(apply(FI[, 1:(mplus - 1), drop = FALSE], 2, sd)^2)
  std_W <- var_T / (mn_var_FI * (mplus - 2))
  W <- rbind(
    matrix(rnorm((mplus - 1) * D), mplus - 1, D) %*%
      diag(c(sqrt(std_W)), D),
    matrix(0, 1, D)
  )
  W[mplus, ] <- colMeans(TT) - colMeans(FI %*% W)
  W
}

#' Beta Initialization from Nearest-Neighbor Distances
#'
#' @param Y Numeric matrix (K x D), mapped points.
#' @return Numeric scalar, initial beta.
#' @keywords internal
.gtm_bi <- function(Y) {
  y_inter_dist <- .gtm_dist(Y, Y) + diag(1e10, nrow(Y))
  mean_nn <- mean(apply(y_inter_dist, 2, min))
  2 / mean_nn
}

#' Sort Responsibilities
#'
#' @param R Numeric matrix (K x N), responsibilities.
#' @return Sorted responsibility matrix.
#' @keywords internal
.gtm_sort <- function(R) {
  idx <- colSums(R^2)
  R[, order(idx)]
}

#' E-Step: Compute Responsibilities and Log-Likelihood
#'
#' Computes posterior responsibilities R and log-likelihood given
#' squared distance matrix and inverse variance beta.
#'
#' @param DIST Numeric matrix (K x N), squared distances.
#' @param minDist Numeric vector (N), column minima of DIST.
#' @param maxDist Numeric vector (N), column maxima of DIST.
#' @param beta Numeric > 0, inverse variance.
#' @param D Integer, data dimensionality.
#' @param md Integer, computation mode (0, 1, or 2).
#' @return List with llh (scalar log-likelihood) and R (K x N matrix).
#' @keywords internal
.gtm_resp <- function(DIST, minDist, maxDist, beta, D, md) {
  if (md < 0 || md > 2 || md != floor(md)) {
    stop("Unknown mode of calculation")
  }
  if (is.matrix(beta) || is.matrix(D)) {
    stop("beta and D should be scalars - mismatch of arguments?")
  }
  if (D != floor(D)) {
    stop("Invalid value for D")
  }
  K <- nrow(DIST)
  N <- ncol(DIST)
  if (md > 0) {
    dist_corr <- (maxDist + minDist) / 2
    dist_corr <- pmin(dist_corr, minDist + 700 * (2 / beta))
    # Bug fix: vectorized with sweep instead of for loop
    DIST <- sweep(DIST, 2, dist_corr, "-")
  }
  R <- exp((-beta / 2) * DIST)
  if (md < 2) {
    r_sum <- colSums(R)
  } else {
    r_sum <- colSums(.gtm_sort(R))
  }
  # Bug fix: vectorized with sweep instead of for loop
  R <- sweep(R, 2, r_sum, "/")
  if (md < 1) {
    llh <- sum(log(r_sum)) + N * ((D / 2) * log(beta / (2 * pi)) - log(K))
  } else {
    llh <- sum(log(r_sum) + dist_corr * (-beta / 2)) +
      N * ((D / 2) * log(beta / (2 * pi)) - log(K))
  }
  list(llh = llh, R = R)
}

#' Wrapper for .gtm_resp in Mode 0
#'
#' @param DIST Numeric matrix (K x N).
#' @param beta Numeric > 0.
#' @param D Integer.
#' @return List with llh and R.
#' @keywords internal
.gtm_resp3 <- function(DIST, beta, D) {
  .gtm_resp(DIST, beta, D, beta, D, 0)
}

#' Setup 1D Latent Space
#'
#' Creates latent grid, basis functions, and initializes weights for a
#' 1D GTM.
#'
#' @param TT Numeric matrix (N x D), training data.
#' @param numsamp Integer, number of latent sample points.
#' @param numbasis Integer, number of basis function centers.
#' @param s Numeric > 0, basis width multiplier.
#' @return List with X, MU, FI, W, beta, sigma.
#' @keywords internal
.gtm_stp1 <- function(TT, numsamp, numbasis, s) {
  if (floor(numsamp) != numsamp || floor(numbasis) != numbasis ||
      numsamp < 0 || numbasis < 0) {
    stop("Incorrect arguments")
  }
  if (s <= 0) {
    stop("Argument s must have strict positive value")
  }
  X <- .gtm_pts(numsamp)
  MU <- .gtm_pts(numbasis)
  MU <- MU * (numbasis / (numbasis - 1))
  sigma <- s * (MU[2] - MU[1])
  FI <- .gtm_gbf(MU, sigma, X)
  pci_result <- .gtm_pci_beta(TT, X, FI)
  list(X = X, MU = MU, FI = FI, W = pci_result$W, beta = pci_result$beta,
       sigma = sigma)
}

#' Setup 2D Latent Space
#'
#' Creates 2D latent grid (rectangular or hexagonal), basis functions,
#' and initializes weights.
#'
#' @param TT Numeric matrix (N x D), training data.
#' @param numsamp Grid dimensions for sample points. Scalar (square grid) or
#'   length-2 vector (xdim, ydim).
#' @param numbasis Grid dimensions for basis centers. Scalar (square grid) or
#'   length-2 vector.
#' @param s Numeric > 0, basis width multiplier.
#' @param kind "rect" or "hex".
#' @return List with X, MU, FI, W, beta, sigma.
#' @keywords internal
.gtm_stp2 <- function(TT, numsamp, numbasis, s, kind = "rect") {
  if (s <= 0) {
    stop("Argument s must have strict positive value")
  }
  if (length(numsamp) == 1 && length(numbasis) == 1) {
    xdim <- sqrt(numsamp)
    ydim <- xdim
    fxdim <- sqrt(numbasis)
    fydim <- fxdim
    if (xdim != floor(xdim) || fxdim != floor(fxdim)) {
      stop("Invalid number of basis functions or latent variable size.")
    }
  } else {
    xdim <- numsamp[1]
    ydim <- numsamp[2]
    # BUG FIX: was numsamp[1]/numsamp[2], should be numbasis[1]/numbasis[2]
    fxdim <- numbasis[1]
    fydim <- numbasis[2]
  }
  if (kind == "rect") {
    X <- .gtm_rctg(xdim, ydim)
    MU <- .gtm_rctg(fxdim, fydim)
  } else {
    X <- .gtm_hex(xdim, ydim)
    MU <- .gtm_hex(fxdim, fydim)
  }
  MU <- MU * (fxdim / (fxdim - 1))
  sigma <- s * abs(MU[1, 1] - MU[2, 1])
  FI <- .gtm_gbf(MU, sigma, X)
  pci_result <- .gtm_pci_beta(TT, X, FI)
  list(X = X, MU = MU, FI = FI, W = pci_result$W, beta = pci_result$beta,
       sigma = sigma)
}

#' EM Training for GTM
#'
#' Runs the EM algorithm for the specified number of cycles.
#'
#' @param TT Numeric matrix (N x D), training data.
#' @param FI Numeric matrix (K x (M+1)), basis function activations.
#' @param W Numeric matrix ((M+1) x D), initial weight matrix.
#' @param l Numeric >= 0, regularization parameter.
#' @param cycles Integer, number of EM iterations.
#' @param beta Numeric > 0, initial inverse variance.
#' @param md Integer, distance mode (default 1).
#' @param quiet Logical, suppress printing (default FALSE).
#' @param min_sing Numeric, minimum singular value threshold (default 0.01).
#' @param tol Numeric, convergence tolerance on log-likelihood change
#'   (default 0, meaning no early stopping).
#' @return List with W, beta, loglik (vector of log-likelihoods per cycle),
#'   R (final responsibilities), converged (logical).
#' @importFrom stats sd
#' @keywords internal
.gtm_trn <- function(TT, FI, W, l, cycles, beta, md = 1, quiet = FALSE,
                     min_sing = 0.01, tol = 0) {
  loglik <- numeric(cycles)
  FI_T <- t(FI)
  K <- nrow(FI)
  mplus <- ncol(FI)
  N <- nrow(TT)
  D <- ncol(TT)
  ND <- N * D
  A <- matrix(0, mplus, mplus)
  if (l > 0) {
    LAMBDA <- l * rep(1, mplus)
    LAMBDA[mplus] <- 0
  }
  gtm_dist_result <- .gtm_dist(TT, FI %*% W, md)
  if (md > 0) {
    gtm_min_dist <- gtm_dist_result$minDist
    gtm_max_dist <- gtm_dist_result$maxDist
    gtm_global_dist <- gtm_dist_result$DIST
  } else {
    gtm_global_dist <- gtm_dist_result
  }
  R_final <- NULL
  converged <- FALSE
  n_iter <- cycles

  for (cycle in 1:cycles) {
    resp_result <- .gtm_resp(gtm_global_dist, gtm_min_dist, gtm_max_dist,
                             beta, D, md)
    R_final <- resp_result$R
    llh <- resp_result$llh
    loglik[cycle] <- llh
    if (!quiet) {
      cat(sprintf("Cycle: %d\tlogLH: %g\tBeta: %g\n", cycle, llh, beta))
    }
    # Check convergence
    if (tol > 0 && cycle > 1) {
      if (abs(loglik[cycle] - loglik[cycle - 1]) < tol) {
        converged <- TRUE
        n_iter <- cycle
        loglik <- loglik[1:cycle]
        break
      }
    }
    # M-step
    if (l > 0) {
      A <- t(t(FI_T) * rowSums(R_final)) %*% FI +
        diag(LAMBDA / beta, nrow = mplus)
    } else {
      A <- t(t(FI_T) * rowSums(R_final)) %*% FI
    }
    chol_result <- chol(A, pivot = TRUE)
    if (attr(chol_result, "rank") != nrow(A)) {
      # Bug fix: use warning() instead of print()
      if (!quiet) {
        warning("M-Step matrix singular, using pseudo-inverse.")
      }
      svd_result <- svd(A)
      sv <- svd_result$d
      # Bug fix: use n_valid_sv instead of N to avoid shadowing
      n_valid_sv <- sum(sv > min_sing)
      if (!quiet) {
        warning(sprintf("Using %d out of %d eigenvalues",
                        n_valid_sv, nrow(A)))
      }
      if (n_valid_sv < 1) {
        stop("Very singular matrix in M-step")
      }
      svd_inverse <- svd_result$v[, 1:n_valid_sv, drop = FALSE] %*%
        diag(1 / sv[1:n_valid_sv], nrow = n_valid_sv) %*%
        t(svd_result$u[, 1:n_valid_sv, drop = FALSE])
      W <- svd_inverse %*% (FI_T %*% (R_final %*% TT))
    } else {
      oo <- order(attr(chol_result, "pivot"))
      W <- chol2inv(chol_result)[oo, oo] %*% (FI_T %*% (R_final %*% TT))
    }
    # Update distances
    gtm_dist_result <- .gtm_dist(TT, FI %*% W, md)
    if (md > 0) {
      gtm_min_dist <- gtm_dist_result$minDist
      gtm_max_dist <- gtm_dist_result$maxDist
      gtm_global_dist <- gtm_dist_result$DIST
    } else {
      gtm_global_dist <- gtm_dist_result
    }
    # Bug fix: update beta from trained values, not just initial
    beta <- ND / sum(colSums(gtm_global_dist * R_final))
  }
  list(W = W, beta = beta, loglik = loglik, R = R_final,
       converged = converged, n_iter = n_iter)
}
