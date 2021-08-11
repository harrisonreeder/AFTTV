#' Generate Matrix of Values of Integrated Piecewise Constant Function
#'
#' This helper function takes in a vector of event times, and
#'   generates a matrix having each row showing the time accrued within each
#'   interval between consecutive elements of a vector of knots.
#'
#' @param y vector of event times.
#' @param knots increasing vector of cutpoints. If it does not start with 0, one will be appended to the start.
#'   However, it should not include Inf at the end.
#' @param intercept if true, includes column corresponding with 'first' interval, and if false does not.
#'   It makes sense to include an intercept if the time-varying covariate is not also included in the "baseline",
#'   otherwise, there would be an identifiability issue.
#'
#' @return a numeric matrix with, with (number of nonzero knots + 1) columns, and with rows corresponding to elements of y.
#' @export
pw_cum_mat <- function(y, knots, intercept=TRUE){
  if(knots[1] != 0){knots <- c(0,knots)}
  #vector giving the length of the intervals between each knot
  knots_diff <- diff(knots)
  #count of the number of intervals in each list
  num_int <- c(length(knots))
  n <- length(y)
  #matrix with each row being an observation, and each column being an interval, capturing how much time is accrued in each interval
  knots_mat <- matrix(c(knots_diff,0),nrow=n,ncol=num_int,byrow = TRUE)
  #an n length vector of integers indicating which interval the observation is in
  cut_cats <- findInterval(x = y, vec = knots)
  #an n length vector capturing the residual time accrued in the final interval of each observation
  y_last <- y-knots[cut_cats]
  #a loop through each observation to finalize the matrix of time intervals
  for(i in 1:n){
    knots_mat[i,cut_cats[i]] <- y_last[i]
    knots_mat[i,-c(1:cut_cats[i])] <- 0
  }

  #removing the intercept changes the parameterization, so that every parameter is a change from the "reference"
  #group which is the first interval. (it basically subtracts the parameter from the first interval from every subsequent interval
  #could be worth changing to instead be the 'chained' parameterization,
  #in which we have each one be difference from the one previous?
  if(!intercept){
    if(ncol(knots_mat) <= 2){stop("include at least two nonzero knots.")}
    knots_mat <- knots_mat[,-1]
  }
  return(knots_mat)
}
