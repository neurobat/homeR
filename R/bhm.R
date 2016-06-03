#' Bayesian Heating Model
#' 
#' Estimates the parameters of a building's heating model.
#' 
#' \code{bhm} assumes that the heating energy for a building has been measured
#' over several time periods (not necessarily of equal length). The \code{data}
#' data frame should have one row per measurement period. The energy vector (whose
#' name is given on the left-hand side of the formula) will have the total energy
#' measured during each period. The daily temperature vector (whose name is given on 
#' the right-hand side of the formula) will have either a vector of average daily
#' temperatures (when each measurement period is just one day) or a list of vectors
#' (when each measurement period can be an arbitrary number of days).
#' 
#' @param formula an object of class "\link{formula}": a description of which variable
#' holds the energy readouts and which variable holds the daily temperatures.
#' @param data a data frame in which the energy and daily temperatures are to be found.
#' @return \code{bhm} returns an object of \link{class} "\code{bhm}". The generic
#' accessor functions \code{coefficients} and \code{residuals} extract the usual
#' information from the fitted model.
#' @examples 
#' set.seed(1111)
#' 
#' # Simple, but unrealistic parameters
#' K <- 1
#' tb <- 1
#' DHW <- 1
#' sigma <- 1e-2
#' temps <- tb + c(-2, -1, 0, 1)
#' 
#' # With daily measurements
#' E <- K * pmax(tb - temps, 0) + DHW + rnorm(length(temps), 0, sigma)
#' fourDayData <- data.frame(E = E, T = temps)
#' fourDayData
#' fit <- bhm(E ~ T, fourDayData)
#' coef(fit)
#' resid(fit)
#' 
#' # With two-day measurements
#' fourTimesTwoDayData <- with(fourDayData,
#'                             data.frame(E = 2 * E,
#'                             T = I(lapply(T, function(x) c(x, x)))))
#' fit2 <- bhm(E ~ T, fourTimesTwoDayData)
#' coef(fit2)
#' resid(fit2)
#' @export
bhm <- function(formula, data) {
  Q <- energy(formula, data)
  stopifnot(all(Q >= 0))
  T <- temperatures(formula, data)
  fit <- posteriorMode(Q, T, minusLogPosterior)
  result <- list()
  result$coefficients <- fit@coef
  result$residuals <- with(as.list(fit@coef),
                           Q - epochEnergy(K, tb, DHW, T))
  class(result) <- "bhm"
  result
}

energy <- function(formula, data) {
  terms <- terms(formula)
  eval(attr(terms, "variables"), envir = data)[[attr(terms, "response")]]
}

temperatures <- function(formula, data) {
  terms <- terms(formula)
  data[[attr(terms, "term.labels")]]
}

posteriorMode <- function(Q, T, logd) {
  fit <- stats4::mle(logd(Q, T), method = "Nelder")
  # restart once
  fit <- stats4::mle(logd(Q, T), method = "Nelder",
                     start = as.list(fit@coef), control = list(maxit = 1000))
  fit
}

minusLogPosterior <- function(energy, temperatures) {
  function(K = 1, tb = 1, DHW = 1, sigma = 1) {
    if (sigma <= 0 || tb < 0 || K <= 0)
      Inf
    else {
      predictedEnergy <- epochEnergy(K, tb, DHW, temperatures)
      stopifnot(length(predictedEnergy) == length(energy))
      epochLength <- sapply(temperatures, length)
      (1 + length(energy)) * log(sigma) + 
        chisquare(energy, predictedEnergy, sqrt(epochLength) * sigma) / 2
    }
  }
}

chisquare <- function(energy, predictedEnergy, epochSigma) {
  stopifnot(length(energy) == length(predictedEnergy))
  stopifnot(length(energy) == length(epochSigma))
  sum(((energy - predictedEnergy) / epochSigma) ^ 2)
}

epochEnergy <- function(K, tb, DHW, temperatures) {
  result <- sapply(temperatures, function(temperature) sum(dailyEnergy(K, tb, DHW, temperature)))
  stopifnot(length(result) == length(temperatures))
  result
}

dailyEnergy <- function(K, tb, DHW, temperature) K * pos(tb - temperature) + DHW

pos <- function(x) pmax(x, 0)

residuals.bhm <- function(bhm) bhm$residuals