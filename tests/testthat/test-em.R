test_that("1D EM training converges with increasing log-likelihood", {
  TT <- matrix(3:61 / 20, ncol = 1)
  TT <- cbind(TT, TT + 1.25 * sin(2 * TT))
  setup <- gtm:::.gtm_stp1(TT, 20, 5, 2)
  result <- gtm:::.gtm_trn(TT, setup$FI, setup$W, l = 0, cycles = 20,
                            beta = setup$beta, quiet = TRUE)
  # Log-likelihood should generally increase (monotone for EM)
  ll <- result$loglik
  # Check that most steps increase (allow tiny numerical dips)
  increases <- diff(ll) >= -1e-8
  expect_true(sum(increases) >= length(increases) - 1)
  # Beta should be positive

  expect_true(result$beta > 0)
  # W should have correct dimensions
  expect_equal(nrow(result$W), ncol(setup$FI))
  expect_equal(ncol(result$W), ncol(TT))
})

test_that("2D EM training works (key bug fix)", {
  TT <- matrix(3:61 / 20, ncol = 1)
  TT <- cbind(TT, TT + 1.25 * sin(2 * TT))
  setup <- gtm:::.gtm_stp2(TT, 81, 25, 2)
  result <- gtm:::.gtm_trn(TT, setup$FI, setup$W, l = 0, cycles = 20,
                            beta = setup$beta, quiet = TRUE)
  expect_true(result$beta > 0)
  expect_equal(nrow(result$W), ncol(setup$FI))
  expect_equal(ncol(result$W), ncol(TT))
  # Log-likelihood should be finite
  expect_true(all(is.finite(result$loglik)))
})

test_that("EM training with regularization works", {
  TT <- matrix(3:61 / 20, ncol = 1)
  TT <- cbind(TT, TT + 1.25 * sin(2 * TT))
  setup <- gtm:::.gtm_stp1(TT, 20, 5, 2)
  result <- gtm:::.gtm_trn(TT, setup$FI, setup$W, l = 0.1, cycles = 10,
                            beta = setup$beta, quiet = TRUE)
  expect_true(result$beta > 0)
  expect_true(all(is.finite(result$loglik)))
})

test_that("EM training with convergence tolerance stops early", {
  TT <- matrix(3:61 / 20, ncol = 1)
  TT <- cbind(TT, TT + 1.25 * sin(2 * TT))
  setup <- gtm:::.gtm_stp1(TT, 20, 5, 2)
  result <- gtm:::.gtm_trn(TT, setup$FI, setup$W, l = 0, cycles = 100,
                            beta = setup$beta, quiet = TRUE, tol = 0.01)
  # Should converge before 100 iterations
  expect_true(result$n_iter < 100)
  expect_true(result$converged)
})

test_that("responsibilities sum to 1 per column", {
  TT <- matrix(3:61 / 20, ncol = 1)
  TT <- cbind(TT, TT + 1.25 * sin(2 * TT))
  setup <- gtm:::.gtm_stp1(TT, 20, 5, 2)
  result <- gtm:::.gtm_trn(TT, setup$FI, setup$W, l = 0, cycles = 5,
                            beta = setup$beta, quiet = TRUE)
  col_sums <- colSums(result$R)
  expect_equal(col_sums, rep(1, ncol(result$R)), tolerance = 1e-10)
})
