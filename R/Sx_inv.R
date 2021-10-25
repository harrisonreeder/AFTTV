#' Compute Inverse Covariate-Adjusted Survival Curve
#'
#' @param para
#' @param p_obj
#' @param x_tv
#' @param beta_tv
#' @param x_base
#' @param beta_base
#' @param baseline
#' @param tv_type
#' @param inv_basis
#' @param knots
#' @param deg
#' @param tstar
#' @param accelfactor
#' @param logSinv
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
Sx_inv <- function(para, p_obj, x_tv=NULL, beta_tv=NULL,
                   x_base=NULL, beta_base=NULL, baseline="weibull",
                   tv_type, inv_basis=NULL, knots=NULL,
                   accelfactor=TRUE, logSinv = FALSE){

  # browser()

  #assume that x_base is non-null
  nP_base <- ncol(x_base)
  nP0 <- if(tolower(baseline)=="weibull") 2 else 0 #for now, hardcode weibull

  #FOR NOW, HARDCODE WEIBULL DISTRIBUTION BUT NOT FOREVER!
  if(tolower(baseline)=="weibull"){
    intercept_temp <- para[1]
    logsigma_temp <- para[2]
    S0_inv <- function(x){(-log(x))^(exp(-logsigma_temp)) * exp(intercept_temp)}
  } else{ stop("only weibull model for now...")}

  S0_p <- S0_inv(p_obj)

  beta_base_temp <- if(nP_base > 0) para[(1+nP0):(nP0+nP_base)] else 0

  #for now, just assume any additional terms are time-varying betas
  beta_tv_temp <- if(tolower(tv_type) !="baseline") para[-(1:(nP0+nP_base))] else 0

  #compute inverse V
  Vinv_temp <- Vx_inv(t_obj = S0_p,
                      beta_base = beta_base_temp, x_base = x_base,
                      beta_tv = beta_tv_temp, x_tv = x_tv,
                      tv_type=tv_type, knots=knots)

  if(accelfactor){Vinv_temp <- Vinv_temp/S0_p}
  if(logSinv){Vinv_temp <- log(Vinv_temp)}
  Vinv_temp

}
