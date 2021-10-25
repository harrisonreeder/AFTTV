#' Compute Integrated Covariate Process V
#'
#' This function computes the integrated covariate process V used for time-varying AFT models.
#'
#' @param t_obj
#' @param x_base
#' @param beta_base
#' @param x_tv
#' @param beta_tv
#' @param xbeta_base
#' @param xbeta_tv
#' @param tv_type
#' @param basis
#' @param knots
#' @param ...
#'
#' @return
#' @export
Vx <- function(t_obj, x_base=NULL, beta_base=NULL, x_tv=NULL, beta_tv=NULL,
               xbeta_base=NULL, xbeta_tv=NULL,
               tv_type, basis=NULL, knots=NULL){
  #browser()

  #First, work through what to do with the time-varying inputs
  #If given xbeta_tv, use it!
  #otherwise, compute it from x_tv and beta_tv
  #if given beta_tv only, set x_tv = 1 (for plotting purposes)
  if(is.null(xbeta_tv)){
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
  }

  stopifnot(tv_type %in% c("baseline", "constant", "step", "logplusone", "piecewise", "other"))
  temp <- NULL
  if(tv_type=="baseline"){
    temp <- t_obj
  } else if(tv_type=="constant"){
    temp <- t_obj*exp(-xbeta_tv)
  } else if(tv_type=="step"){
    temp <- t_obj + (exp(-xbeta_tv) - 1) * (t_obj-tstar) * as.numeric(t_obj>tstar)
  } else if(tv_type=="logplusone"){

    flag <- xbeta_tv != 1

    temp <- rep(NA, length(flag))
    temp[flag] <- ((t_obj+1)^(1-xbeta_tv)-1)/(1-xbeta_tv)
    temp[!flag] <- log(1+t_obj)

    #use fast "ifelse" type statement (https://stackoverflow.com/questions/38004924/speeding-up-ifelse-without-writing-c-c)
    # flag <- xbeta_tv != 1
    # temp <- flag * ((t_obj+1)^(1-xbeta_tv)-1)/(1-xbeta_tv) + (!flag) * log(1+t_obj)
  } else if(tv_type=="piecewise"){
    #need just the "diagonal" of the matrix product of the time basis, and the exp(-x beta\trans) matrix
    #https://stackoverflow.com/questions/42569698/how-to-just-calculate-the-diagonal-of-a-matrix-product-in-r

    #this formula parameterizes beta_tv's as three different (for use with 'intercept' basis)
    # temp <- as.numeric(rowSums(t_obj * exp(-xbeta_tv)))

    if(is.null(basis)){
      if(is.null(knots)){stop("Need to supply a basis, or a vector of knots.")}
      basis <- pw_cum_mat(y=t_obj, knots=knots,intercept = FALSE)
    }

    #this formula instead parameterizes "relative" to the first interval
    #so, first beta_tv is effect in second interval relative to first, second is in third interval relative to first, etc.
    temp <- as.numeric(rowSums(basis * (exp(-xbeta_tv)-1)) + t_obj)

  } else if(tv_type=="other"){
    stopifnot(!is.null(basis))
    #figure out how to do numerical integration of some kind here
  }

  #Lastly, work out what to do with the basesline covariates
  if(is.null(xbeta_base)){
    if(is.null(x_base)){
      if(is.null(beta_base)){
        warning("we're assuming no baseline covariates, because xbeta_base, x_base, and beta_base are all NULL.")
        xbeta_base <- 0
      } else{
        warning("x_base and xbeta_base are both NULL, so we assume all baseline covariates are set to 1.")
        xbeta_base <- sum(beta_base)
      }
    } else{
      #calculate xbeta_base using x_base and beta_base
      xbeta_base <- x_base %*% as.matrix(beta_base)
    }
  }

  as.numeric(temp * exp(-xbeta_base))
}
