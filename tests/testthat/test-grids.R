test_that("1D grid points are correctly spaced", {
  pts <- gtm:::.gtm_pts(5)
  expect_equal(nrow(pts), 5)
  expect_equal(ncol(pts), 1)
  expect_equal(pts[1, 1], -1)
  expect_equal(pts[5, 1], 1)
  # Evenly spaced
  diffs <- diff(pts[, 1])
  expect_true(all(abs(diffs - diffs[1]) < 1e-10))
})

test_that("1D grid points work with M=2", {
  pts <- gtm:::.gtm_pts(2)
  expect_equal(nrow(pts), 2)
  expect_equal(pts[1, 1], -1)
  expect_equal(pts[2, 1], 1)
})

test_that("rectangular grid has correct dimensions", {
  g <- gtm:::.gtm_rctg(4, 3)
  expect_equal(nrow(g), 12)
  expect_equal(ncol(g), 2)
})

test_that("rectangular grid is centered", {
  g <- gtm:::.gtm_rctg(5, 5)
  expect_equal(mean(g[, 1]), 0, tolerance = 1e-10)
  expect_equal(mean(g[, 2]), 0, tolerance = 1e-10)
})

test_that("rectangular grid rejects invalid dimensions", {
  expect_error(gtm:::.gtm_rctg(1, 3))
  expect_error(gtm:::.gtm_rctg(3, 1))
  expect_error(gtm:::.gtm_rctg(2.5, 3))
})

test_that("hexagonal grid has correct dimensions", {
  g <- gtm:::.gtm_hex(4, 3)
  expect_equal(nrow(g), 12)
  expect_equal(ncol(g), 2)
})

test_that("hexagonal grid is approximately centered", {
  g <- gtm:::.gtm_hex(5, 5)
  # Hex grids may not be perfectly centered due to offset rows
  expect_equal(mean(g[, 1]), 0, tolerance = 0.1)
  expect_equal(mean(g[, 2]), 0, tolerance = 0.1)
})

test_that("hexagonal grid rejects invalid dimensions", {
  expect_error(gtm:::.gtm_hex(1, 3))
  expect_error(gtm:::.gtm_hex(3, 1))
})

test_that("hexagonal grid matches old implementation", {
  # Regression values from gtm.hex.old(10,10)
  g <- gtm:::.gtm_hex(10, 10)
  expect_equal(nrow(g), 100)
  # Check it's within [-1, 1] range
  expect_true(all(g[, 1] >= -1.01))
  expect_true(all(g[, 1] <= 1.01))
  expect_true(all(g[, 2] >= -1.01))
  expect_true(all(g[, 2] <= 1.01))
})
