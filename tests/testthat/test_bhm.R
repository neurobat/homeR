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

expect_within <- function(object, expected, pm)
  expect_equal(as.vector(object), expected, tolerance = pm, scale = 1)

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

test_that("works with synthetic data from paper", {
  load("fakeMonthlyEnergy.RData")
  model <- bhm(Energy ~ DailyMeans, fakeMonthlyEnergy)
  coefs <- coef(model)
  stderr <- sqrt(diag(vcov(model)))
  expect_within(coefs['K'], 20, pm = 5e-2)
  expect_within(stderr['K'], 0.1, pm = 0.1)
  expect_within(coefs['tb'], 12, pm = 5e-2)
  expect_within(stderr['tb'], 0.2, pm = 0.02)
  expect_within(coefs['DHW'], 100, pm = 2)
  expect_within(stderr['DHW'], 2, pm = 1)
  expect_within(coefs['sigma'], 30, pm = 1)
  expect_within(stderr['sigma'], 3, pm = 1)
})

test_that("works with heat counter data", {
  load("alv.RData") # heat counter data from paper
  model <- bhm(Energy ~ DailyMeans, alv)
  coefs <- coef(model)
  expect_within(coefs['K'], 2.4, pm = 0.1)
  expect_within(coefs['tb'], 17.7, pm = 0.1)
  expect_within(coefs['DHW'], 0.7, pm = 0.1)
  expect_within(coefs['sigma'], 7.8, pm = 0.1)
})

test_that("provided log-posterior has maximum at estimated coefficients", {
  model <- bhm(E ~ T, fourDayData)
  coefs <- coef(model)
  nonOptimalCoefs <- lapply(coefs, function(x) jitter(rep(x, 2)))
  logp <- logposterior(model)
  logPosteriorMode <- do.call(logp, as.list(coefs))
  logPosteriorElsewhere <- do.call(logp, nonOptimalCoefs)
  expect_true(all(logPosteriorElsewhere < logPosteriorMode))
  # with a small jitter, the logp should never decrease by more than, say, 20
  expect_true(all(logPosteriorElsewhere > logPosteriorMode - 20))
})

test_that("user can specify an arbitrary DHW", {
  threeDaysBelowBaseTemp <- fourDayData[1:3, ]
  model <- bhm(E ~ T, data = threeDaysBelowBaseTemp, baseLoad = 0)
  coefs <- coef(model)
  expect_within(coefs['K'], K, pm = 0.1)
  expect_within(coefs['DHW'], 0, pm = 0.1)
})