library(homeR)

library(plyr)
context("Bayesian heating model")

# The heating energy model assumes that the heating energy grows linearly as
# the outdoor temperature decreases, but remains constant above the base temperature.
# For unit tests we use a model with base temp = 1, K = 1, DHW = 1 and sigma = 1e-2.
# 
#        ^
#        |
#        |               |
#        |      \
#        |       \
#        |        \      |
# Energy |         \
#        |          \
#        |           \   |
#        |            \
#        | total heat loss
#        |    coef (K)  \|
#        |             ( \---------------
#        +--------------------------------->
#                        |          Outdoor
#                      base          temp.
#                      temp.

set.seed(1111)

K <- 1
tb <- 1
DHW <- 1
sigma <- 1e-2

temps <- tb + c(-2, -1, 0, 1)
E <- K * pmax(tb - temps, 0) + DHW + rnorm(length(temps), 0, sigma)

fourDayData <- data.frame(E = E, T = temps)

expect_equal_coefs <- function(model)
  expect_equal(coef(model), c(K = K, tb = tb, DHW = DHW, sigma = sigma), tolerance = 1e-2, scale = 1)

test_that("finds coefs on four-day data", {
  model <- bhm(E ~ T, fourDayData)
  expect_equal_coefs(model)
  expect_lt(sum(resid(model)^2), 1e-1)
})

test_that("negative energy yields an error", {
  dataWithNegativeEnergy <- transform(fourDayData, E = -E)
  expect_error(bhm(E ~ T, dataWithNegativeEnergy))
})

test_that("provides residuals", {
  model <- bhm(E ~ T, fourDayData)
  expect_length(residuals(model), nrow(fourDayData))
})

test_that("find coefs on four times two-day data", {
  fourTwoDayEpochData <- adply(fourDayData, 1,
                               summarise, E = 2 * E, T = I(list(c(T, T))))
  model <- bhm(E ~ T, fourTwoDayEpochData)
  expect_equal_coefs(model)
})

test_that("works on epochs of different legths", {
  load("oilConsumption.RData") # load the `goodPeriods' data frame
  model <- bhm(Energy ~ DailyMeans, goodPeriods)
  coefs <- coef(model)
  expect_equal(coefs['K'], c(K = 6.1), tolerance = 5e-2)
  expect_equal(coefs['tb'], c(tb = 17.7), tolerance = 5e-2)
  expect_equal(coefs['DHW'], c(DHW = 0), tolerance = 2, scale = 1)
  expect_equal(coefs['sigma'], c(sigma = 100), tolerance = 5, scale = 1)
})