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