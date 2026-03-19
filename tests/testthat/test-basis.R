test_that("Gaussian basis has correct dimensions", {
  MU <- gtm:::.gtm_pts(4)
  X <- gtm:::.gtm_pts(10)
  FI <- gtm:::.gtm_gbf(MU, sigma = 1, X)
  expect_equal(nrow(FI), 10)  # K
  expect_equal(ncol(FI), 5)   # M + 1 (4 basis + 1 bias)
})

test_that("Gaussian basis has bias column of ones", {
  MU <- gtm:::.gtm_pts(3)
  X <- gtm:::.gtm_pts(8)
  FI <- gtm:::.gtm_gbf(MU, sigma = 1, X)
  expect_true(all(FI[, ncol(FI)] == 1))
})

test_that("Gaussian basis values are in (0, 1]", {
  MU <- gtm:::.gtm_pts(5)
  X <- gtm:::.gtm_pts(20)
  FI <- gtm:::.gtm_gbf(MU, sigma = 0.5, X)
  # Excluding bias column
  FI_body <- FI[, -ncol(FI)]
  expect_true(all(FI_body > 0))
  expect_true(all(FI_body <= 1))
})

test_that("Gaussian basis rejects dimension mismatch", {
  MU <- matrix(rnorm(6), ncol = 2)  # 3x2
  X <- matrix(rnorm(15), ncol = 3)  # 5x3
  expect_error(gtm:::.gtm_gbf(MU, sigma = 1, X))
})

test_that("linear basis has correct dimensions", {
  X <- gtm:::.gtm_pts(10)
  FI <- gtm:::.gtm_lbf(X)
  expect_equal(nrow(FI), 10)
  expect_equal(ncol(FI), 2)  # L + 1
  expect_true(all(FI[, 2] == 1))
})

test_that("2D Gaussian basis works", {
  MU <- gtm:::.gtm_rctg(3, 3)  # 9x2
  X <- gtm:::.gtm_rctg(5, 5)   # 25x2
  FI <- gtm:::.gtm_gbf(MU, sigma = 1, X)
  expect_equal(nrow(FI), 25)
  expect_equal(ncol(FI), 10)  # 9 + 1
})
