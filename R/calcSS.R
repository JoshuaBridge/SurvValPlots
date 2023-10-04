#' Wrapper for pmsampsize with r-squared estimation
#'
#' @import pmsampsize
#'
#' @param rate Event rate at time
#' @param time Time point to show performance at
#' @param parameters Number of parameters in the model
#' @param meanfollowup The mean follow-up time
#' @param imprecision Allowance for imprecision in the rate and mean follow-up
#'
#' @return GGPlot object
#' @export
#'
#' @examples
#' calcSS(
#'   rate = 0.1,
#'   time = 3,
#'   parameters = 5,
#'   meanfollowup = 4
#' )
calcSS = function(
    rate,
    time,
    parameters,
    meanfollowup,
    imprecision=0.2
){

  rate = rate-rate*imprecision
  meanfollowup = meanfollowup-meanfollowup*imprecision
  lnL = rate*meanfollowup*(log(rate*meanfollowup)-1)
  maxR2 = 1-exp((2*lnL))
  estR2 = maxR2*0.15
  pmsampsize("s",
             parameters=parameters,
             rsquared=estR2,
             rate=rate,
             timepoint=time,
             meanfup=meanfollowup)
}

