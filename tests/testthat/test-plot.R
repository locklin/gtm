# Smoke tests: ensure all plot types run without error

test_that("plot convergence works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "convergence"))
  dev.off()
})

test_that("plot means works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "means"))
  dev.off()
})

test_that("plot means with labels works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "means", labels = iris$Species))
  dev.off()
})

test_that("plot modes works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "modes", labels = iris$Species))
  dev.off()
})

test_that("plot density works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "density"))
  dev.off()
})

test_that("plot umatrix works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "umatrix"))
  dev.off()
})

test_that("plot magnification works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "magnification"))
  dev.off()
})

test_that("plot components single variable works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "components", component = 1))
  dev.off()
})

test_that("plot components all variables works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "components"))
  dev.off()
})

test_that("plot hex grid density works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), grid_type = "hex", n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "density"))
  dev.off()
})

test_that("plot 1D model convergence works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), latent_dim = 1, n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "convergence"))
  dev.off()
})

test_that("plot 1D model density works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), latent_dim = 1, n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "density"))
  dev.off()
})

test_that("plot with property coloring works", {
  data(iris)
  m <- GtmFit(as.matrix(iris[, 1:4]), n_iter = 10)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(m, type = "means", property = iris$Sepal.Length))
  dev.off()
})
