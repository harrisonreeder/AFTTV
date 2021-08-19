#' Simulate Data under Time-Varying AFT Model
#'
#' @param x_base
#' @param x_tv
#' @param beta_base_true
#' @param beta_tv_true
#' @param tv_type
#' @param S0_inv_func
#' @param cens
#' @param knots
#' @param deg
#' @param tstar
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sim_AFTtv <- function (x_base, x_tv, beta_base_true, beta_tv_true,
                        tv_type, S0_inv_func, cens, knots, deg=1,tstar=1e100,...) {
  n <- dim(x_base)[1]
  p <- dim(x_base)[2]

  #this is just one way of sampling from the baseline distribution,
  #though you could also sample from it "directly" !!
  u <- stats::runif(n = n, min = 0, max = 1)
  T_baseline <- S0_inv_func(u, ...)

  T_temp <- Vx_inv(t_obj = T_baseline, x_base=x_base, x_tv=x_tv,
                   beta_base=beta_base_true, beta_tv=beta_tv_true,
                   tv_type=tv_type, knots=knots, deg=deg,tstar=tstar,...)

  delta <- rep(NA, n)
  y <- T_temp
  C_temp <- stats::runif(n = n, min = cens[1], max = cens[2])

  ind1 <- which(T_temp < C_temp)
  y[ind1] <- T_temp[ind1]
  delta[ind1] <- 1
  ind0 <- which(T_temp >= C_temp)
  y[ind0] <- C_temp[ind0]
  delta[ind0] <- 0

  data.frame(y=y, delta=delta)
}




#' Simulate Data under Time-Varying AFT Model From Select Distributions
#'
#' @param x_base
#' @param x_tv
#' @param beta_base_true
#' @param beta_tv_true
#' @param tv_type
#' @param dist
#' @param intercept
#' @param scale
#' @param cens
#' @param knots
#' @param deg
#' @param tstar
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sim_AFTtv_dist <- function (x_base, x_tv, beta_base_true, beta_tv_true,
                       tv_type, dist, intercept, scale, cens, knots, deg=1,tstar=1e100,...) {
  n <- dim(x_base)[1]
  p <- dim(x_base)[2]

  #this is just one way of sampling from the baseline distribution,
  #though you could also sample from it "directly" !!
  if(dist == "weibull"){
    T_baseline <- rweibull(n = n,scale=exp(intercept),shape=1/scale)
  } else if(dist == "lognormal"){
    T_baseline <- rlnorm(n = n,meanlog = intercept,sdlog=scale)
  }

  T_temp <- Vx_inv(t_obj = T_baseline, x_base=x_base, x_tv=x_tv,
                   beta_base=beta_base_true, beta_tv=beta_tv_true,
                   tv_type=tv_type, knots=knots, deg=deg,tstar=tstar,...)

  delta <- rep(NA, n)
  y <- T_temp
  C_temp <- stats::runif(n = n, min = cens[1], max = cens[2])

  ind1 <- which(T_temp < C_temp)
  y[ind1] <- T_temp[ind1]
  delta[ind1] <- 1
  ind0 <- which(T_temp >= C_temp)
  y[ind0] <- C_temp[ind0]
  delta[ind0] <- 0

  data.frame(y=y, delta=delta)
}



