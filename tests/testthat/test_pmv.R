library(homeR)
context("Predicted Mean Vote")

test_that("pmv() handles scalars", {
  vote <- pmv(clo = 1, met = 1.2, air.temp = 19, saturation = 40)
  expect_gt(vote, -0.6)
  expect_lt(vote, -0.5)
})

test_that("pmv() handles vectors", {
  votes <- pmv(clo = 1, met = 1.2, air.temp = c(19, 30), sat = 40)
  expect_length(votes, 2)
  expect_equal(votes, c(pmv(clo = 1, met = 1.2, air.temp = 19, sat = 40), pmv(clo = 1, 
    met = 1.2, air.temp = 30, sat = 40)))
}) 
