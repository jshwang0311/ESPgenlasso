context("Exact Solution Path algorithm for genlasso")
library(GenLassoTest)

ex.y <- matrix(c(1000), nrow = 1)
ex.X <- matrix(c(1,1,0),nrow = 1)
ex.D <- matrix(c(1,0,0,0,1,1,0,0,-1), nrow = 3)
object <- ESPgenlasso(ex.y,ex.X,ex.D,genlasso.option = FALSE)
test_that("Testing ESPgenlasso", {
  expect_equal(length(object), 11)
  expect_equal(sum(abs(object$beta[,1])), 0)
})

ex.y <- matrix(c(1000), nrow = 1)
ex.X <- matrix(c(1,1,0),nrow = 1)
ex.D <- matrix(c(1,0,0,0,1,1,0,0,-1), nrow = 3)
object <- ESPgenlasso(ex.y,ex.X,ex.D,genlasso.option = FALSE)
test_that("Testing get.beta", {
  expect_equal(sum(abs(get.beta(object, 2000))), 0)
  expect_equal(get.beta(object, 700)[1], 200)
  expect_equal(get.beta(object, 700)[2], 100)
  expect_equal(get.beta(object, 700)[3], 100)
})

ex.y <- matrix(c(1000), nrow = 1)
ex.X <- matrix(c(1,1,0),nrow = 1)
ex.D <- matrix(c(1,0,0,0,1,1,0,0,-1), nrow = 3)
object <- ESPgenlasso(ex.y,ex.X,ex.D,genlasso.option = FALSE)
test_that("Testing get.u", {
  expect_equal(get.u(object, 700)[1], 700)
  expect_equal(get.u(object, 700)[2], 700)
  expect_equal(get.u(object, 700)[3], 0)
})
