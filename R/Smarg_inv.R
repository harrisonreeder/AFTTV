#' Compute inverse of regression-standardized survival function
#'
#' @param p
#' @param x_base_aug
#' @param x_tv_aug
#' @param beta_base
#' @param beta_tv
#' @param int
#' @param shp
#' @param tv_type
#' @param knots
#' @param lower
#' @param upper
#'
#' @return
#' @export
Smarg_inv = function(p,x_base_aug, x_tv_aug, beta_base, beta_tv,
                     int, shp, tv_type, knots=NULL,lower=0,upper=100){
  # browser()
  tryCatch(stats::uniroot(f = function (t){p -
      Smarg(t=t, x_base_aug = x_base_aug, x_tv_aug=x_tv_aug,
            beta_base = beta_base, beta_tv = beta_tv, int = int, shp = shp,
            tv_type = tv_type, knots = knots)},
      lower = lower, upper = upper)$root,
      error=function(e){return(NA)})
}
