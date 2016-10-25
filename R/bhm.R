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
#' @param baseLoad a optional constant base load, e.g. for domestic hot water preparation.
#' @return \code{bhm} returns an object of \link{class} "\code{bhm}". The generic
#' accessor functions \code{coefficients}, \code{vcov} and \code{residuals} extract the usual
#' information from the fitted model, while \code{logposterior} will return a function
#' that evaluates the log-posterior as a function of the parameters.
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
#' \dontrun{
#' fit <- bhm(E ~ T, fourDayData)
#' coef(fit)
#' resid(fit)
#' }
#' 
#' # With two-day measurements
#' fourTimesTwoDayData <- with(fourDayData,
#'                             data.frame(E = 2 * E,
#'                             T = I(lapply(T, function(x) c(x, x)))))
#' fit2 <- bhm(E ~ T, fourTimesTwoDayData)
#' coef(fit2)
#' resid(fit2)
#' @export
bhm <- function(formula, data, baseLoad = NULL) {
  Q <- energy(formula, data)
  stopifnot(all(Q >= 0))
  T <- temperatures(formula, data)
  result <- list()
  result$logPosterior <- logPosteriorFunction(Q, T)
  fit <- posteriorMode(result$logPosterior, baseLoad)
  stopifnot('DHW' %in% names(fit$par) || !is.null(baseLoad))
  result$coefficients <- c(fit$par, DHW = baseLoad)
  result$residuals <- with(as.list(result$coefficients),
                           Q - epochEnergy(K, tb, DHW, T))
  result$vcov <- solve(fit$hessian)
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

posteriorMode <- function(logposterior, baseLoad = NULL) {
  obj <- function(par, ...) {
   - do.call(logposterior, as.list(c(par, ...)))
  }
  start <- list(K = 1, tb = 10, DHW = 1, sigma = 1)
  fixed <- if (!is.null(baseLoad)) list(DHW = baseLoad) else list()
  start[names(fixed)] <- NULL
  fit <- stats::optim(start, obj, gr = NULL, fixed, method = "SANN", control = list(maxit = 20000))
  fit <- stats::optim(fit$par, obj, gr = NULL, fixed, control = list(maxit = 1000), hessian = TRUE)
  fit
}

logPosteriorFunction <- function(energy, temperatures) {
  mlp <- function(K = 10, tb = 20, DHW = 0, sigma = 100) {
    if (sigma <= 0) return(-Inf)
    predictedEnergy <- epochEnergy(K, tb, DHW, temperatures)
    stopifnot(length(predictedEnergy) == length(energy))
    epochLength <- sapply(temperatures, length)
    - ((1 + length(energy)) * log(sigma) + 
         chisquare(energy, predictedEnergy, sqrt(epochLength) * sigma) / 2)
  }
  Vectorize(mlp)
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

#' @export
residuals.bhm <- function(object, ...) object$residuals

#' @export
vcov.bhm <- function(object, ...) object$vcov

#' Log-posterior of a Bayesian Heating Model
#' 
#' Provides the log-posterior of a heating model given the data, as a function
#' of the model's parameters.
#' 
#' @param bhm a fitted model returned by a call to \code{bhm()}
#' @return a function of the model's parameters (currently K, tb, DHW and sigma)
#' @export
logposterior <- function(bhm) bhm$logPosterior