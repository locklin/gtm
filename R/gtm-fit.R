#' Fit a Generative Topographic Map
#'
#' Fits a GTM model to high-dimensional data using the EM algorithm
#' from Bishop, Svensen, and Williams (1998). The GTM is a probabilistic
#' alternative to Self-Organizing Maps (SOM) that provides a principled
#' density model in data space via a smooth nonlinear mapping from a
#' low-dimensional latent space.
#'
#' @param data Numeric matrix (N x D), the high-dimensional data.
#' @param latent_dim Integer, 1 or 2 (default 2).
#' @param grid_size Grid dimensions. Scalar (square grid) or length-2 vector
#'   for 2D. Default: `ceiling(sqrt(5 * sqrt(nrow(data))))`.
#' @param n_basis Number of basis functions. Scalar or length-2 vector.
#'   Default: `ceiling(grid_size / 2)`.
#' @param basis_width Numeric > 0. Width of Gaussian basis functions relative
#'   to inter-center distance. Default: 2.
#' @param grid_type `"rect"` or `"hex"` (2D only). Default: `"rect"`.
#' @param n_iter Integer, maximum number of EM iterations. Default: 100.
#' @param lambda Numeric >= 0, regularization. Default: 0.01.
#' @param init `"pca"` or `"random"`. Default: `"pca"`.
#' @param tol Numeric, convergence tolerance on log-likelihood. Default: 1e-6.
#' @param verbose Logical. Default: FALSE.
#'
#' @return S3 object of class `"gtm"` with components:
#' \describe{
#'   \item{W}{Weight matrix ((M+1) x D).}
#'   \item{beta}{Scalar inverse variance.}
#'   \item{FI}{Basis function activations (K x (M+1)).}
#'   \item{X}{Latent grid points (K x L).}
#'   \item{MU}{Basis function centers (M x L).}
#'   \item{sigma}{Basis function width.}
#'   \item{Y}{Mapped points in data space (K x D).}
#'   \item{R}{Final responsibility matrix (K x N).}
#'   \item{loglik}{Vector of log-likelihoods per iteration.}
#'   \item{data}{Training data matrix (N x D).}
#'   \item{latent_dim}{Latent dimensionality (1 or 2).}
#'   \item{grid_size}{Grid dimensions used.}
#'   \item{grid_type}{"rect" or "hex".}
#'   \item{n_iter}{Number of iterations actually run.}
#'   \item{converged}{Logical, did training converge?}
#' }
#'
#' @references
#' Bishop, C. M., Svensen, M., & Williams, C. K. I. (1998).
#' GTM: The Generative Topographic Mapping.
#' *Neural Computation*, 10(1), 215-234.
#' \doi{10.1162/089976698300017953}
#'
#' @examples
#' # Fit a 2D GTM to the iris data
#' data(iris)
#' model <- GtmFit(as.matrix(iris[, 1:4]), verbose = TRUE, n_iter = 20)
#' print(model)
#'
#' # 1D latent space
#' model_1d <- GtmFit(as.matrix(iris[, 1:4]), latent_dim = 1, n_iter = 20)
#'
#' @export
GtmFit <- function(data, latent_dim = 2L, grid_size = NULL,
                   n_basis = NULL, basis_width = 2, grid_type = "rect",
                   n_iter = 100L, lambda = 0.01, init = "pca",
                   tol = 1e-6, verbose = FALSE) {
  # Validate inputs
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("'data' must be a numeric matrix.")
  }
  N <- nrow(data)
  D <- ncol(data)
  if (N < 2) stop("Need at least 2 data points.")
  if (D < 1) stop("Data must have at least 1 column.")
  latent_dim <- as.integer(latent_dim)
  if (!(latent_dim %in% 1:2)) {
    stop("'latent_dim' must be 1 or 2.")
  }
  grid_type <- match.arg(grid_type, c("rect", "hex"))
  init <- match.arg(init, c("pca", "random"))
  if (basis_width <= 0) stop("'basis_width' must be positive.")
  if (n_iter < 1) stop("'n_iter' must be at least 1.")
  if (lambda < 0) stop("'lambda' must be non-negative.")

  # Default grid size
  if (is.null(grid_size)) {
    gs <- ceiling(sqrt(5 * sqrt(N)))
    if (gs < 2) gs <- 2
    if (latent_dim == 1L) {
      grid_size <- gs^2  # total number of points for 1D
    } else {
      grid_size <- c(gs, gs)
    }
  } else {
    if (latent_dim == 1L) {
      grid_size <- as.integer(grid_size[1])
      if (grid_size < 2) stop("'grid_size' must be >= 2.")
    } else {
      if (length(grid_size) == 1) {
        grid_size <- c(as.integer(grid_size), as.integer(grid_size))
      } else {
        grid_size <- as.integer(grid_size[1:2])
      }
      if (any(grid_size < 2)) stop("'grid_size' dimensions must be >= 2.")
    }
  }

  # Default n_basis
  if (is.null(n_basis)) {
    if (latent_dim == 1L) {
      n_basis <- max(2L, ceiling(grid_size / 2))
    } else {
      n_basis <- pmax(2L, ceiling(grid_size / 2))
    }
  } else {
    if (latent_dim == 1L) {
      n_basis <- as.integer(n_basis[1])
    } else {
      if (length(n_basis) == 1) {
        n_basis <- c(as.integer(n_basis), as.integer(n_basis))
      } else {
        n_basis <- as.integer(n_basis[1:2])
      }
    }
  }

  # Setup latent space and basis functions
  if (latent_dim == 1L) {
    setup <- .gtm_stp1(data, grid_size, n_basis, basis_width)
  } else {
    setup <- .gtm_stp2(data, grid_size, n_basis, basis_width,
                       kind = grid_type)
  }

  X <- setup$X
  MU <- setup$MU
  FI <- setup$FI
  sigma <- setup$sigma

  if (init == "pca") {
    W <- setup$W
    beta <- setup$beta
  } else {
    W <- .gtm_ri(data, FI)
    beta <- .gtm_bi(FI %*% W)
  }

  # Train
  trained <- .gtm_trn(data, FI, W, l = lambda, cycles = n_iter,
                      beta = beta, md = 1, quiet = !verbose,
                      tol = tol)

  # Build result object
  Y <- FI %*% trained$W

  structure(list(
    W         = trained$W,
    beta      = trained$beta,
    FI        = FI,
    X         = X,
    MU        = MU,
    sigma     = sigma,
    Y         = Y,
    R         = trained$R,
    loglik    = trained$loglik,
    data      = data,
    latent_dim = latent_dim,
    grid_size = grid_size,
    grid_type = grid_type,
    n_iter    = trained$n_iter,
    converged = trained$converged
  ), class = "gtm")
}
