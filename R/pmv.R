#' Predicted Mean Vote
#' 
#' Computes Fanger's predicted mean vote
#' 
#' Compute the predicted mean vote for one or more combinations of 
#' clo, met, air temperature and moisture saturation. The inputs arguments 
#' can be scalars or vectors.
#' 
#' @param clo Thermal insulation of clothing in [clo] (underwear, blouse/shirt, slacks/trousers, jacket, socks and shoes are approximately 1 clo)
#' @param met Physical activity in [met] (one person seated at rest is approximately 1 met)
#' @param air.temp Indoor air temperature (assumed equal to mean radiant temperature) in [C]
#' @param saturation Ratio of moisture content to moisture content of saturated air at the same temperature, in [\%] (approximately the same thing as relative humidity)
#' @return The predicted mean vote, a value between -3 (cold) to +3 (hot)
#' @references CIBSE Guide A, section 1.4 and 1.A1.2 (from which this implementation is derived)
#' @examples 
#' # With scalars
#' pmv(clo=1.0,
#'     met=1.2,
#'     air.temp=19,
#'     saturation=40)
#' # With vectors
#' pmv(clo=c(1.0, 1.5),
#'     met=c(1.2, 0.6),
#'     air.temp=c(19, 30),
#'     sat=c(35, 40))
#' @export
pmv <- function(clo, met, air.temp, saturation) {
  # Adapted from CIBSE Guide A, 1.A1.2, p. 51
  saturated.vapour.pressure <- exp(16.6536-4030.183/(air.temp+235)) # kPa
  vapour.pressure <- saturated.vapour.pressure*10*saturation # Pa
  clothing.insulation <- 0.155*clo
  met.rate <- met*58.15
  external.work <- 0 # simplification
  int.heat.prod <- met.rate - external.work
  clo.area.fact <- ifelse(clothing.insulation<0.078,
                          1 + 1.29*clothing.insulation,
                          1.05 + 0.645*clothing.insulation)
  convective.heat <- 0 # assumes null air velocity
  air.temp.kel <- air.temp + 273
  mean.rad.temp <- air.temp.kel # simplification
  # compute clothing surface temperature by iteration
  clo.surf.temp <- air.temp.kel +
    (35.5 - air.temp)/(3.5*(6.45*clothing.insulation + 0.1))
  p1 <- clothing.insulation*clo.area.fact
  p2 <- p1*3.96
  p3 <- p1*100 # Error in CIBSE Guide
  p4 <- p1*air.temp.kel
  p5 <- 308.7 - 0.028*int.heat.prod + p2*(mean.rad.temp/100)^4
  xn <- clo.surf.temp/100
  xf <- xn
  eps <- 0.000015
  hc <- 0
  repeat {
    xf <- (xf + xn)/2
    heat.transf.nat <- 2.38*abs(100*xf-air.temp.kel)^0.25
    hc <- pmax(convective.heat, heat.transf.nat)
    xn <- (p5 + p4*hc - p2*xf^4)/(100 + p3*hc)
    if (all(abs(xn-xf)<eps)) break
  }
  clo.surf.temp <- 100*xn - 273
  # end iteration
  heat.loss.skin <- 3.05*0.001*(5733-6.99*int.heat.prod-vapour.pressure) # CIBSE Guide has an error here
  heat.loss.sweat <- pmax(0.42*(int.heat.prod - 58.15), 0)
  heat.loss.respiration <- 1.7*0.00001*met.rate*(5867 - vapour.pressure)
  heat.loss.dry <- 0.0014*met.rate*(34 - air.temp)
  heat.loss.rad <- 3.96*clo.area.fact*(xn^4 - (mean.rad.temp/100)^4)
  heat.loss.conv <- clo.area.fact*hc*(clo.surf.temp - air.temp)
  therm.sens <- 0.303*exp(-0.036*met.rate) + 0.028
  therm.sens*(int.heat.prod - heat.loss.skin - heat.loss.sweat
              - heat.loss.respiration - heat.loss.dry
              - heat.loss.rad - heat.loss.conv)
}
pmv <- Vectorize(pmv)
