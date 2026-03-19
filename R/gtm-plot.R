#' Plot Method for GTM
#'
#' Produces a variety of visualizations for a fitted Generative Topographic
#' Map, following the `kohonen::plot.kohonen` pattern of `type=` dispatch.
#'
#' @param x A fitted `"gtm"` object from [GtmFit()].
#' @param type Type of plot. One of:
#' \describe{
#'   \item{`"means"`}{Data points projected to R-weighted mean latent positions.}
#'   \item{`"modes"`}{Data points projected to highest-responsibility node.}
#'   \item{`"density"`}{Heatmap of cumulative responsibility per node.}
#'   \item{`"umatrix"`}{Heatmap of neighbor distances between mapped points.}
#'   \item{`"components"`}{Per-variable heatmap of R-weighted data values.}
#'   \item{`"magnification"`}{Heatmap of local area distortion (Jacobian).}
#'   \item{`"convergence"`}{Line plot of log-likelihood over iterations.}
#' }
#' @param labels Optional factor or character vector (length N) for coloring
#'   points in `"means"` and `"modes"` plots.
#' @param property Optional numeric vector (length N) for coloring points
#'   by a continuous variable in `"means"` and `"modes"` plots.
#' @param component Integer or character, which variable(s) to show for
#'   `"components"` type. Default shows all (one per panel).
#' @param main Plot title. If NULL, a default is used.
#' @param palette Color palette: a function(n) returning n colors, or a
#'   character vector of colors. Default: `hcl.colors(64, "YlOrRd")`.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @examples
#' data(iris)
#' model <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 20)
#' plot(model, type = "convergence")
#'
#' @importFrom graphics plot plot.new plot.window points axis box title
#'   polygon rect par mtext lines legend
#' @importFrom grDevices hcl.colors
#' @export
plot.gtm <- function(x, type = c("means", "modes", "density", "umatrix",
                                  "components", "magnification", "convergence"),
                     labels = NULL, property = NULL, component = NULL,
                     main = NULL, palette = NULL, ...) {
  type <- match.arg(type)
  model <- x

  switch(type,
    means = .plot_projection(model, "means", labels, property, main,
                             palette, ...),
    modes = .plot_projection(model, "modes", labels, property, main,
                             palette, ...),
    density = .plot_density(model, main, palette, ...),
    umatrix = .plot_umatrix(model, main, palette, ...),
    components = .plot_components(model, component, main, palette, ...),
    magnification = .plot_magnification(model, main, palette, ...),
    convergence = .plot_convergence(model, main, ...)
  )
  invisible(model)
}

#' @keywords internal
.plot_projection <- function(model, proj_type, labels, property, main,
                             palette, ...) {
  proj <- predict(model, type = proj_type)
  if (is.null(main)) {
    main <- paste0("GTM ", proj_type, " projection")
  }

  L <- model$latent_dim
  if (L == 2L) {
    # Determine colors
    if (!is.null(labels)) {
      labels <- as.factor(labels)
      n_levels <- nlevels(labels)
      if (is.null(palette)) {
        pal <- grDevices::hcl.colors(n_levels, "Set 2")
      } else if (is.function(palette)) {
        pal <- palette(n_levels)
      } else {
        pal <- palette
      }
      cols <- pal[as.integer(labels)]
    } else if (!is.null(property)) {
      if (is.null(palette)) {
        pal <- grDevices::hcl.colors(64, "Viridis")
      } else if (is.function(palette)) {
        pal <- palette(64)
      } else {
        pal <- palette
      }
      pr <- range(property, na.rm = TRUE)
      if (pr[1] == pr[2]) {
        col_idx <- rep(32, length(property))
      } else {
        col_idx <- floor((property - pr[1]) / (pr[2] - pr[1]) *
                           (length(pal) - 1)) + 1
        col_idx <- pmin(pmax(col_idx, 1), length(pal))
      }
      cols <- pal[col_idx]
    } else {
      cols <- "steelblue"
    }

    graphics::plot(proj[, 1], proj[, 2], col = cols, pch = 16, cex = 0.8,
                   xlab = "Latent dim 1", ylab = "Latent dim 2",
                   main = main, asp = 1, ...)
    if (!is.null(labels)) {
      graphics::legend("topright", legend = levels(labels),
                       col = pal[seq_len(nlevels(labels))],
                       pch = 16, cex = 0.7, bg = "white")
    }
  } else {
    # 1D
    if (!is.null(labels)) {
      labels <- as.factor(labels)
      n_levels <- nlevels(labels)
      if (is.null(palette)) {
        pal <- grDevices::hcl.colors(n_levels, "Set 2")
      } else if (is.function(palette)) {
        pal <- palette(n_levels)
      } else {
        pal <- palette
      }
      cols <- pal[as.integer(labels)]
    } else {
      cols <- "steelblue"
    }
    graphics::stripchart(proj[, 1] ~ seq_len(nrow(proj)), method = "jitter",
                         col = cols, pch = 16,
                         xlab = "Latent position", ylab = "",
                         main = main, ...)
  }
}

#' @keywords internal
.plot_density <- function(model, main, palette, ...) {
  if (is.null(main)) main <- "GTM node density"
  rho <- rowSums(model$R)
  .draw_grid_heatmap(model, rho, main = main, palette = palette)
}

#' @keywords internal
.plot_umatrix <- function(model, main, palette, ...) {
  if (is.null(main)) main <- "GTM U-matrix"
  u_vals <- .gtm_umatrix(model)
  .draw_grid_heatmap(model, u_vals, main = main, palette = palette)
}

#' @keywords internal
.plot_components <- function(model, component, main, palette, ...) {
  D <- ncol(model$data)
  R <- model$R
  r_sums <- rowSums(R)
  # Weighted mean of each data variable per node
  # comp_j[k] = sum_n(R[k,n] * data[n,j]) / sum_n(R[k,n])

  if (is.null(component)) {
    component <- seq_len(D)
  } else if (is.character(component)) {
    component <- match(component, colnames(model$data))
    if (any(is.na(component))) {
      stop("'component' names not found in data columns.")
    }
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  n_comp <- length(component)
  if (n_comp > 1) {
    nc <- ceiling(sqrt(n_comp))
    nr <- ceiling(n_comp / nc)
    graphics::par(mfrow = c(nr, nc), mar = c(2, 2, 3, 1))
  }

  cnames <- colnames(model$data)
  for (j in component) {
    comp_vals <- (R %*% model$data[, j]) / r_sums
    if (is.null(main)) {
      ttl <- if (!is.null(cnames)) cnames[j] else paste("Variable", j)
    } else {
      ttl <- main
    }
    .draw_grid_heatmap(model, comp_vals, main = ttl, palette = palette)
  }
}

#' @keywords internal
.plot_magnification <- function(model, main, palette, ...) {
  if (is.null(main)) main <- "GTM magnification factors"
  mag <- .gtm_magnification(model)
  .draw_grid_heatmap(model, mag, main = main, palette = palette)
}

#' @keywords internal
.plot_convergence <- function(model, main, ...) {
  if (is.null(main)) main <- "GTM training convergence"
  n <- length(model$loglik)
  graphics::plot(seq_len(n), model$loglik, type = "b", pch = 16, cex = 0.5,
                 xlab = "Iteration", ylab = "Log-likelihood",
                 main = main, ...)
}
