#' Positive Part
#' 
#' Computes the positive part of each element of a vector.
#' 
#' @param x A numeric vector
#' @return The argument with all negative entries replaced by 0
pos <- function(x) pmax(x, 0)

#' Daily Energy Model
#' 
#' Computes the daily energy model without the noise component.
#' 
#' @param dailyMeans A vector of mean daily temperatures
#' @param K The building's total heat loss coefficient
#' @param baseTemp The building's base temperature
#' @param DHW The building's average base load
#' @return A vector of predicted energy consumption for each day.
dailyEnergyModel <- function(dailyMeans, K = 20, baseTemp = 12, DHW = 100)
  K * pos(baseTemp - dailyMeans) + DHW

#' Period Energy Model
#' 
#' Computes the energy consumption over several days, without the noise component.
#' 
#' @inheritParams dailyEnergyModel
#' @return The total predicted energy consumption for that period.
periodEnergyModel <- function(...) sum(dailyEnergyModel(...))

#' Chi-square for a single datum
#' 
#' Computes the chi-square statistic for a single reporting period.
#' 
#' @param response The value predicted by the model given the data and the parameters
#' @param model The model that predicts the response given data and parameters
#' @param sigma The standard deviation that should be used for computing the Chi-square
#' @param ... Extra arguments passed to the model
#' @return The chi-square for that single datum.
chisquare1 <- function(response, model, sigma, ...) {
  ((model(...) - response)/sigma)^2
}

#' Chi-square for several reporting periods
#' 
#' @param periods A data frame with two elements: `Energy`, a vector of the energy
#' recorded during each period, and `DailyMeans`, a list of mean daily temperature
#' vectors for each reporting period.
#' @param U, baseTemp, DHW, sigma As above
chisquare <- function(periods, U, baseTemp, DHW, sigma) {
  stopifnot("Energy" %in% names(periods))
  stopifnot("DailyMeans" %in% names(periods))
  chisquares <- apply(periods, MARGIN = 1, function(period) chisquare1(period$Energy, 
    periodEnergyModel, sqrt(length(period$DailyMeans)) * sigma, 
    period$DailyMeans, U, baseTemp, DHW))
  sum(chisquares)
}

#' The log-likelihood of the data given the parameters
#' 
#' @param periods, U, baseTemp, DHW, sigma As above
loglikelihood <- function(periods, U, baseTemp, DHW, sigma) {
  if (sigma <= 0) 
    return(-Inf)
  if (DHW < 0) 
    return(-Inf)
  -(1 + nrow(periods)) * log(sigma) - chisquare(periods, U, baseTemp, 
    DHW, sigma)/2
}
loglikelihood <- Vectorize(loglikelihood, c("U", "baseTemp", "DHW", 
  "sigma"))

#' Estimate building parameters and their uncertainties from the data
#' 
#' @param periods The energy and mean outdoor temperature, as above
#' @return A list of two elements: `params`, a list of estimated parameters, 
#' and `covariance`, the covariance matrix
estimateParameters <- function(periods) {
  # Find parameters that maximise log-likelihood
  objective <- function(params, data) {
    curried <- function(...) -loglikelihood(data, ...)
    do.call(curried, as.list(params))
  }
  optimal <- optim(par = c(10, 20, 0, 100), fn = objective, gr = NULL, 
    periods, hessian = TRUE)  # for covariance matrix
  params = list(U = 0, baseTemp = 0, DHW = 0, sigma = 0)
  params[] <- optimal$par
  list(params = params, covariance = solve(optimal$hessian))
} 
