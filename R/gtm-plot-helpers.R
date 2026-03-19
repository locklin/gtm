# Internal helpers for plot.gtm

#' Draw a Filled Hexagon
#'
#' @param cx,cy Center coordinates.
#' @param r Radius (circumradius).
#' @param col Fill color.
#' @param border Border color.
#' @keywords internal
.draw_hexagon <- function(cx, cy, r, col, border = NA) {
  angles <- seq(0, 2 * pi, length.out = 7)
  xs <- cx + r * cos(angles)
  ys <- cy + r * sin(angles)
  graphics::polygon(xs, ys, col = col, border = border)
}

#' Draw Grid Heatmap (Rect or Hex)
#'
#' Maps a vector of values to colors and draws them on the GTM latent grid.
#'
#' @param model A gtm object.
#' @param values Numeric vector (length K), one value per grid node.
#' @param main Plot title.
#' @param palette Color palette function or character vector.
#' @keywords internal
.draw_grid_heatmap <- function(model, values, main = "", palette = NULL) {
  X <- model$X
  K <- nrow(X)

  if (is.null(palette)) {
    pal <- grDevices::hcl.colors(64, "YlOrRd", rev = TRUE)
  } else if (is.function(palette)) {
    pal <- palette(64)
  } else {
    pal <- palette
  }

  # Map values to colors
  val_range <- range(values, na.rm = TRUE)
  if (val_range[1] == val_range[2]) {
    col_idx <- rep(32, K)
  } else {
    col_idx <- floor((values - val_range[1]) /
                       (val_range[2] - val_range[1]) * (length(pal) - 1)) + 1
    col_idx <- pmin(pmax(col_idx, 1), length(pal))
  }
  cols <- pal[col_idx]

  if (model$latent_dim == 2L) {
    # Determine cell radius
    dists <- as.matrix(stats::dist(X))
    diag(dists) <- Inf
    min_dist <- min(dists)
    r <- min_dist * 0.55

    x_range <- range(X[, 1])
    y_range <- range(X[, 2])
    pad <- r * 1.5

    graphics::plot.new()
    graphics::plot.window(xlim = x_range + c(-pad, pad),
                          ylim = y_range + c(-pad, pad),
                          asp = 1)

    if (model$grid_type == "hex") {
      for (i in seq_len(K)) {
        .draw_hexagon(X[i, 1], X[i, 2], r, col = cols[i], border = "grey80")
      }
    } else {
      for (i in seq_len(K)) {
        graphics::rect(X[i, 1] - r, X[i, 2] - r,
                       X[i, 1] + r, X[i, 2] + r,
                       col = cols[i], border = "grey80")
      }
    }
    graphics::title(main = main)
    graphics::axis(1)
    graphics::axis(2)
    graphics::box()
  } else {
    # 1D: line plot
    ord <- order(X[, 1])
    graphics::plot(X[ord, 1], values[ord], type = "h",
                   col = cols[ord], lwd = 3,
                   xlab = "Latent position", ylab = "Value",
                   main = main)
  }
}

#' Compute U-Matrix (Neighbor Distances)
#'
#' For each grid node, computes the mean distance to its neighbors
#' in the mapped data space.
#'
#' @param model A gtm object.
#' @return Numeric vector (length K).
#' @keywords internal
.gtm_umatrix <- function(model) {
  Y <- model$Y
  X <- model$X
  K <- nrow(X)

  # Find neighbors: nodes within 1.5 * min_inter_node_distance
  latent_dists <- as.matrix(stats::dist(X))
  diag(latent_dists) <- Inf
  min_dist <- min(latent_dists)
  threshold <- min_dist * 1.5

  mapped_dists <- as.matrix(stats::dist(Y))

  u_vals <- numeric(K)
  for (k in seq_len(K)) {
    neighbors <- which(latent_dists[k, ] < threshold)
    if (length(neighbors) > 0) {
      u_vals[k] <- mean(mapped_dists[k, neighbors])
    }
  }
  u_vals
}

#' Compute Magnification Factors
#'
#' Computes the local area distortion (magnification factor) at each
#' grid node from the Jacobian of the mapping.
#'
#' @param model A gtm object.
#' @return Numeric vector (length K).
#' @keywords internal
.gtm_magnification <- function(model) {
  X <- model$X
  MU <- model$MU
  W <- model$W
  sigma <- model$sigma
  K <- nrow(X)
  L <- ncol(X)
  M <- nrow(MU)
  D <- ncol(W)
  mplus <- nrow(W)  # M+1

  # W_body = W without bias row (M x D)
  W_body <- W[1:M, , drop = FALSE]

  mag <- numeric(K)
  for (k in seq_len(K)) {
    # Compute dPhi/dX at node k: M x L matrix
    # phi_km = exp(-||x_k - mu_m||^2 / (2*sigma^2))
    # dPhi_km/dx_l = -(1/sigma^2) * (x_kl - mu_ml) * phi_km
    diff_km <- sweep(MU, 2, X[k, ], "-")  # M x L: mu_m - x_k
    sq_dist <- rowSums(diff_km^2)
    phi_k <- exp(-sq_dist / (2 * sigma^2))  # M-vector

    # dPhi/dx: M x L, each col l = -(1/sigma^2) * (x_kl - mu_ml) * phi_km
    # = (1/sigma^2) * (mu_ml - x_kl) * phi_km  (since diff = mu - x)
    dPhi <- (1 / sigma^2) * diff_km * phi_k  # M x L

    # Jacobian: D x L = W_body^T %*% dPhi
    J <- t(W_body) %*% dPhi  # D x L

    if (L == 1) {
      mag[k] <- sqrt(sum(J^2))  # ||J||
    } else {
      # sqrt(det(J^T J))
      JtJ <- t(J) %*% J  # L x L
      mag[k] <- sqrt(abs(det(JtJ)))
    }
  }
  mag
}
