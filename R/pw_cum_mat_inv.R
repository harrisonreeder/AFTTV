#' Generate Matrix of Values of Inverse of Integrated Piecewise Constant Function
#'
#' @param y
#' @param knots
#' @param x_base
#' @param beta_base
#' @param x_tv
#' @param beta_tv
#'
#' @return
#' @export
pw_cum_mat_inv <- function(y, knots, x_base=x_base,beta_base=beta_base, x_tv, beta_tv){

  y_adj <- y #as.numeric(y * exp(-x_base %*% beta_base))
  #p0 is number of knots (including 0)

  #these are the knots on the original time scale
  if(knots[1] != 0){knots <- c(0,knots)}

  #n by p0 matrix with knots for each subject, shifted by baseline amount
  adj_knots_mat <- exp(-x_base %*% beta_base) %*% t(knots)
  #n by (p0-1) matrix with differences for each subject, shifted by baseline amount
  # tstar1*exp(-x_base trans beta_base), (tstar2 - tstar1)*exp(-x_base trans beta_base), etc.
  adj_knots_diff <- t(diff(t(adj_knots_mat)))

  #p0
  nP0 <- length(beta_tv)
  n <- length(y)

  stopifnot(nP0 > 1)

  #n by (p0-1) matrix with (exp(-x_i*beta1),exp(-x_i*beta2),...), omitting last beta
  temp_out <- exp(-x_tv %*% t(beta_tv[-nP0]))


  #n by (p0-1) matrix with each column giving corrected interval length between each knot
  #tstar1*exp(-x_base trans beta_base) (tstar2-tstar1)*exp(-x_base trans beta_base - x*beta1)
  adj_knots_diff_mat <- adj_knots_diff
  adj_knots_diff_mat[,-1] <- adj_knots_diff_mat[,-1] * temp_out

  #n by p0 matrix with each column being the jth knot for the ith subject
  adj_knots_mat <- cbind(0,t(apply(adj_knots_diff_mat,MARGIN = 1,cumsum)))

  #add column of 0's as placeholder, as this will be output matrix
  adj_knots_diff_mat <- cbind(adj_knots_diff_mat,0)

  #an n length vector of integers indicating which interval the observation is in
  #https://stackoverflow.com/questions/27472948/vectorization-of-findinterval
  cut_cats <- rowSums(y_adj >= adj_knots_mat) #mild question of whether this should be open or closed

  #an n length vector capturing the residual time accrued in the final interval of each observation
  #https://stackoverflow.com/questions/20617371/getting-elements-of-a-matrix-with-vectors-of-coordinates
  y_last <- y_adj - adj_knots_mat[cbind(1:n,cut_cats)]

  #a loop through each observation to finalize the matrix of time intervals
  for(i in 1:n){
    adj_knots_diff_mat[i,cut_cats[i]] <- y_last[i]
    adj_knots_diff_mat[i,-c(1:cut_cats[i])] <- 0
  }

  #adj_knots_diff_mat tallies up the total accrued time of individual i through all intervals
  #(after transforming to the inverse scale)

  return(adj_knots_diff_mat)
}
