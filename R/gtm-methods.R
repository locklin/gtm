#' Predict Method for GTM
#'
#' Projects data into the latent space of a fitted GTM model.
#'
#' @param object A fitted `"gtm"` object from [GtmFit()].
#' @param newdata Optional numeric matrix (N_new x D). If NULL, uses the
#'   training data stored in the model.
#' @param type One of `"means"` (R-weighted mean latent positions),
#'   `"modes"` (latent position of highest-responsibility node),
#'   or `"responsibilities"` (full K x N responsibility matrix).
#' @param ... Ignored.
#'
#' @return Depends on `type`:
#' \describe{
#'   \item{"means"}{Numeric matrix (N x L) of mean latent positions.}
#'   \item{"modes"}{Numeric matrix (N x L) of mode latent positions.}
#'   \item{"responsibilities"}{Numeric matrix (K x N) of responsibilities,
#'     where K is the number of grid nodes and N is the number of data points.}
#' }
#'
#' @examples
#' data(iris)
#' model <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 20)
#' proj <- predict(model, type = "means")
#' head(proj)
#'
#' @importFrom stats predict
#' @export
predict.gtm <- function(object, newdata = NULL,
                        type = c("means", "modes", "responsibilities"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    R <- object$R
  } else {
    if (!is.matrix(newdata) || !is.numeric(newdata)) {
      stop("'newdata' must be a numeric matrix.")
    }
    if (ncol(newdata) != ncol(object$data)) {
      stop("'newdata' must have the same number of columns as training data.")
    }
    Y <- object$Y
    dist_result <- .gtm_dist(newdata, Y, m = 1)
    resp <- .gtm_resp(dist_result$DIST, dist_result$minDist,
                      dist_result$maxDist, object$beta,
                      ncol(newdata), md = 1)
    R <- resp$R
  }
  X <- object$X
  switch(type,
    means = t(R) %*% X,
    modes = X[apply(R, 2, which.max), , drop = FALSE],
    responsibilities = R
  )
}

#' Print Method for GTM
#'
#' @param x A `"gtm"` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.gtm <- function(x, ...) {
  N <- nrow(x$data)
  D <- ncol(x$data)
  K <- nrow(x$X)
  cat("Generative Topographic Map\n")
  cat(sprintf("  Data: %d x %d\n", N, D))
  cat(sprintf("  Latent dim: %d\n", x$latent_dim))
  if (x$latent_dim == 2L) {
    cat(sprintf("  Grid: %d x %d (%s), %d nodes\n",
                x$grid_size[1], x$grid_size[2], x$grid_type, K))
  } else {
    cat(sprintf("  Grid: %d nodes\n", K))
  }
  cat(sprintf("  Basis functions: %d (+1 bias)\n", ncol(x$FI) - 1))
  cat(sprintf("  Iterations: %d", x$n_iter))
  if (x$converged) cat(" (converged)")
  cat("\n")
  cat(sprintf("  Final log-likelihood: %.4f\n", x$loglik[length(x$loglik)]))
  cat(sprintf("  Final beta: %.4f\n", x$beta))
  invisible(x)
}

#' Summary Method for GTM
#'
#' @param object A `"gtm"` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns a list with summary statistics.
#'
#' @importFrom stats sd
#' @export
summary.gtm <- function(object, ...) {
  N <- nrow(object$data)
  D <- ncol(object$data)
  K <- nrow(object$X)
  means <- predict(object, type = "means")
  out <- list(
    n_data = N,
    n_dims = D,
    n_nodes = K,
    latent_dim = object$latent_dim,
    grid_size = object$grid_size,
    grid_type = object$grid_type,
    n_iter = object$n_iter,
    converged = object$converged,
    final_loglik = object$loglik[length(object$loglik)],
    final_beta = object$beta,
    mean_projection = means
  )
  cat("Generative Topographic Map Summary\n")
  cat("===================================\n")
  cat(sprintf("Data dimensions: %d observations x %d variables\n", N, D))
  cat(sprintf("Latent dimensionality: %d\n", object$latent_dim))
  if (object$latent_dim == 2L) {
    cat(sprintf("Grid: %d x %d %s (%d nodes)\n",
                object$grid_size[1], object$grid_size[2],
                object$grid_type, K))
  } else {
    cat(sprintf("Grid: %d nodes\n", K))
  }
  cat(sprintf("Basis functions: %d (+1 bias)\n", ncol(object$FI) - 1))
  cat(sprintf("Iterations run: %d (converged: %s)\n",
              object$n_iter, object$converged))
  cat(sprintf("Final log-likelihood: %.6f\n",
              object$loglik[length(object$loglik)]))
  cat(sprintf("Final beta (inverse variance): %.6f\n", object$beta))
  # Density stats
  rho <- rowSums(object$R)
  cat(sprintf("Node density: min=%.3f, max=%.3f, sd=%.3f\n",
              min(rho), max(rho), sd(rho)))
  invisible(out)
}
