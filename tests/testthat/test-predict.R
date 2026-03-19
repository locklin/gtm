test_that("predict means returns correct dimensions", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  p <- predict(m, type = "means")
  expect_equal(nrow(p), 150)
  expect_equal(ncol(p), 2)
})

test_that("predict modes returns correct dimensions", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  p <- predict(m, type = "modes")
  expect_equal(nrow(p), 150)
  expect_equal(ncol(p), 2)
})

test_that("predict modes returns valid grid indices", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  modes <- predict(m, type = "modes")
  # All mode positions should be actual grid points
  for (i in seq_len(nrow(modes))) {
    match_found <- any(apply(m$X, 1, function(row) {
      all(abs(row - modes[i, ]) < 1e-10)
    }))
    expect_true(match_found)
  }
})

test_that("predict responsibilities sum to 1", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  R <- predict(m, type = "responsibilities")
  expect_equal(nrow(R), nrow(m$X))
  expect_equal(ncol(R), 150)
  col_sums <- colSums(R)
  expect_equal(col_sums, rep(1, 150), tolerance = 1e-10)
})

test_that("predict with newdata works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  new_data <- as.matrix(iris[1:10, 1:4])
  p <- predict(m, newdata = new_data, type = "means")
  expect_equal(nrow(p), 10)
  expect_equal(ncol(p), 2)
})

test_that("predict with newdata validates dimensions", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  bad_data <- matrix(1:6, ncol = 2)
  expect_error(predict(m, newdata = bad_data))
})

test_that("1D predict returns correct dimensions", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), latent_dim = 1, n_iter = 10)
  p <- predict(m, type = "means")
  expect_equal(nrow(p), 150)
  expect_equal(ncol(p), 1)
})
