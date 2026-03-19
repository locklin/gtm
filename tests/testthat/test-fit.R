test_that("GtmFit returns correct S3 class", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  expect_s3_class(m, "gtm")
})

test_that("GtmFit 2D has all required fields", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  expected_fields <- c("W", "beta", "FI", "X", "MU", "sigma", "Y", "R",
                       "loglik", "data", "latent_dim", "grid_size",
                       "grid_type", "n_iter", "converged")
  for (f in expected_fields) {
    expect_true(f %in% names(m), info = paste("Missing field:", f))
  }
})

test_that("GtmFit 2D dimensions are consistent", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10, grid_size = 6)
  N <- nrow(iris)
  D <- 4
  K <- 36  # 6x6 grid
  M <- ncol(m$FI) - 1
  expect_equal(nrow(m$X), K)
  expect_equal(ncol(m$X), 2)
  expect_equal(nrow(m$FI), K)
  expect_equal(nrow(m$W), M + 1)
  expect_equal(ncol(m$W), D)
  expect_equal(nrow(m$Y), K)
  expect_equal(ncol(m$Y), D)
  expect_equal(nrow(m$R), K)
  expect_equal(ncol(m$R), N)
  expect_equal(nrow(m$data), N)
  expect_equal(ncol(m$data), D)
  expect_equal(m$latent_dim, 2L)
  expect_equal(m$grid_type, "rect")
})

test_that("GtmFit 1D works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), latent_dim = 1, n_iter = 10)
  expect_s3_class(m, "gtm")
  expect_equal(m$latent_dim, 1L)
  expect_equal(ncol(m$X), 1)
})

test_that("GtmFit hex grid works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), grid_type = "hex", n_iter = 10)
  expect_s3_class(m, "gtm")
  expect_equal(m$grid_type, "hex")
})

test_that("GtmFit random init works", {
  set.seed(42)
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), init = "random", n_iter = 10)
  expect_s3_class(m, "gtm")
  expect_true(all(is.finite(m$loglik)))
})

test_that("GtmFit rejects invalid inputs", {
  expect_error(GtmFit(1:10))  # not a matrix
  expect_error(GtmFit(matrix(1:4, 2, 2), latent_dim = 3))
  expect_error(GtmFit(matrix(1, 1, 2)))  # only 1 row
})

test_that("GtmFit convergence detection works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 200, tol = 0.001)
  # With generous tol, should converge well before 200
  if (m$converged) {
    expect_true(m$n_iter < 200)
    expect_equal(length(m$loglik), m$n_iter)
  }
})

test_that("GtmFit with explicit grid_size and n_basis", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), grid_size = c(5, 6),
              n_basis = c(3, 3), n_iter = 5)
  expect_equal(nrow(m$X), 30)  # 5*6
  expect_equal(nrow(m$MU), 9)  # 3*3
})
