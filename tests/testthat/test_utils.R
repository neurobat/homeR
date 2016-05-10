library(homeR)
context("Utils")

model <- list(K = 5, tb = 11, DHW = 10)

test_that("positive part function returns correct values for scalars and vectors", {
  expect_equal(pos(-5), 0)
  expect_equal(pos(1), 1)
  expect_equal(pos(c(-1, 0, 1)), c(0, 0, 1))
})

test_that("daily energy model returns correct values on both sides of base temperature", {
  DHW <- 10; K <- 5; baseTemp <- 11
  model <- dailyEnergyModel(c(-1, 0, 1) + model$tb, K = K, baseTemp = baseTemp, DHW = DHW)
  expect_equal(model, DHW + c(K, 0, 0))
})

test_that("period energy model returns correct value for days on both sides of base temperature", {
  DHW <- 10; K <- 5; baseTemp <- 11
  model <- periodEnergyModel(c(-1, 0, 1) + baseTemp, K = K, baseTemp = baseTemp, DHW = DHW)
  expect_equal(model, K + 3 * DHW)
})

test_that("chi-square1 matches values computed by hand", {
  
})