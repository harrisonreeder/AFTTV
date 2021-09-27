#' Compute regression-standardized survival function
#'
#' @param t
#' @param x_base_aug
#' @param x_tv_aug
#' @param beta_base
#' @param beta_tv
#' @param int
#' @param shp
#' @param tv_type
#' @param knots
#' @param baseline
#' @param ...
#'
#' @return
#' @export
Smarg <- function(t, x_base_aug, x_tv_aug, beta_base, beta_tv,
                  int, shp, tv_type, knots=NULL,baseline="weibull",...){
  #aug prefix is a reminder that the intention is that this is a matrix with
  #the (time-varying) exposure 'set to' a value, and everything else marginalized out
  # browser()
  N <- length(x_tv_aug)
  V_temp <- Vx(t_obj = rep(t,N),
               x_base = x_base_aug, x_tv = x_tv_aug,
               beta_base = beta_base,
               beta_tv = beta_tv,
               tv_type=tv_type,knots=knots,...)
  if(tolower(baseline) %in% c("wb","weibull")){
    S_temp_vec <- stats::pweibull(q=V_temp,
                                  scale = exp(int), shape = shp,lower.tail = FALSE)
  } else{
    S_temp_vec <- stats::plnorm(q=V_temp,
                                  meanlog = int, sdlog = shp,lower.tail = FALSE)
  }
  #assume Weibull for now
  mean(S_temp_vec)
}
