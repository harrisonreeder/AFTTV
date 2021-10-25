

#' Compute Inverse of Integrated Covariate Process V
#'
#' @param t_obj
#' @param x_tv
#' @param beta_tv
#' @param xbeta_base
#' @param x_base
#' @param beta_base
#' @param tv_type
#' @param inv_basis
#' @param knots
#'
#' @return
#' @export
Vx_inv <- function(t_obj, x_tv=NULL, beta_tv=NULL,
                   xbeta_base=NULL, x_base=NULL, beta_base=NULL,
                   tv_type, inv_basis=NULL, knots=NULL){
  # browser()
  temp <- NULL

  if(is.null(xbeta_base)){
    if(is.null(x_base)){
      if(is.null(beta_base)){
        warning("we're assuming no baseline covariates, because xbeta_base, x_base, and beta_base are all NULL.")
        xbeta_base <- 0
        x_base <- numeric(length(t_obj))
        beta_base <- 0
      } else{
        warning("x_base and xbeta_base are both NULL, so we assume all baseline covariates are set to 1.")
        #just assume that we want to plot with x set to 1 for everyone
        xbeta_base <- sum(beta_base)
        x_base <- matrix(1, nrow=length(t_obj), ncol=length(beta_base))
      }
    } else{
      xbeta_base <- x_base %*% as.matrix(beta_base)
    }
  }

  if(tv_type=="piecewise"){
    if(is.null(x_tv)){
      warning("x_tv was null, so we're replacing it with 1's")
      x_tv <- rep(1,length(t_obj))
    }
    if(is.null(inv_basis)){
      inv_basis <- pw_cum_mat_inv(y=t_obj,knots=knots,
                                  x_base=x_base,beta_base=beta_base,
                                  x_tv=x_tv,beta_tv=beta_tv)
    }
    #matrix with columns
    #exp(x_base trans beta_base)  exp(x_base trans beta_base + x_tv * beta_tv1) exp(x_base trans beta_base + x_tv * beta_tv2)
    mult_temp <- cbind(1,exp(x_tv %*% t(beta_tv))) * as.numeric(exp(xbeta_base))

    #vector with row sums of the above quantities
    temp <- as.numeric(rowSums(inv_basis * mult_temp))
  } else{

    #now that we're in non-piecewise
    tvec_xbeta <- exp(xbeta_base) * t_obj

    if(is.null(x_tv)){
      if(is.null(beta_tv)){
        if(tv_type != "baseline"){
          stop("for everything except tv_type='baseline', must provide either xbeta_tv, or beta_tv (and optionally x_tv).")
        }
      } else{
        xbeta_tv <- rep(1,length(t_obj)) %*% t(beta_tv)
      }
    } else{
      xbeta_tv <- x_tv %*% t(beta_tv)
    }


    if(tv_type=="baseline"){temp <- tvec_xbeta}
    else if(tv_type=="constant"){temp <- tvec_xbeta*exp(xbeta_tv)}
    else if(tv_type=="step"){temp <- tvec_xbeta + (exp(xbeta_tv)-1)*(tvec_xbeta-tstar)*as.numeric(tvec_xbeta>tstar)}
    else if(tv_type=="logplusone"){
      #faster ifelse: https://github.com/ICJIA/r-user-group/issues/11
      flag <- xbeta_tv != 1

      temp <- rep(NA, length(flag))
      temp[flag] <- ( (1-xbeta_tv)*tvec_xbeta + 1 )^(1/(1-xbeta_tv)) - 1
      temp[!flag] <- exp(tvec_xbeta) - 1
    } else{stop("tv_type must be 'piecewise', 'baseline', 'constant', 'step', 'logplusone'")}
  }
  as.numeric(temp)
}
