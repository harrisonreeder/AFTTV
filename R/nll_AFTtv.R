#' Title
#'
#' @param para
#' @param y
#' @param delta
#' @param x_base
#' @param x_tv
#' @param baseline
#' @param basis
#' @param tv_type
#' @param deg
#' @param tstar
#' @param ...
#'
#' @return
#' @export
nll_AFTtv <- function (para, y, delta, x_base, x_tv, baseline = "weibull", basis=NULL,
                       tv_type, deg=1,tstar=1e100,...){
  # browser()

  nP_base <- ncol(x_base)
  if(tolower(tv_type)=="baseline"){
    nP_tv <- 0
  } else{
    nP_tv <- if(is.null(basis)) 1 else ncol(basis)
  }
  nP0 <- if(tolower(baseline) %in% c("weibull","lognormal")) 2 else 0
  nP_tot <- nP_base + nP_tv + nP0
  stopifnot(length(para) == nP_tot)

  intercept_temp <- para[1]
  logsigma_temp <- para[2]
  beta_base_temp <- if(nP_base > 0) para[(1+nP0):(nP0+nP_base)] else 0
  beta_tv_temp <- if(nP_tv > 0) para[(1+nP0+nP_base):(nP0+nP_base+nP_tv)] else 0

  xbeta_base_temp <- x_base %*% as.matrix(beta_base_temp)
  if(tolower(tv_type)!="baseline"){
    xbeta_tv_temp <- x_tv %*% t(beta_tv_temp)
  } else{
    xbeta_tv_temp <- rep(0,length(y))
  }

  V_temp <- Vx(t_obj=y, xbeta_base=xbeta_base_temp,
               xbeta_tv=xbeta_tv_temp,basis=basis,
               tv_type=tv_type, deg=deg,tstar=tstar,...)

  if(tolower(baseline)=="weibull"){
    logS0 <- function(x){stats::pweibull(q=x,scale = exp(intercept_temp),
                                  shape = exp(logsigma_temp),
                                  lower.tail = FALSE, log.p = TRUE)}
    logh0 <- function(x){flexsurv::hweibull(x=x,scale = exp(intercept_temp),
                                            shape = exp(logsigma_temp), log = TRUE)}
  } else if(tolower(baseline)=="lognormal"){
    logS0 <- function(x){stats::plnorm(q=x, meanlog = intercept_temp,
                                         sdlog = exp(logsigma_temp),
                                         lower.tail = FALSE, log.p = TRUE)}
    logh0 <- function(x){log(flexsurv::hlnorm(x=x, meanlog = intercept_temp,
                                            sdlog = exp(logsigma_temp)))} #logging manually bc of a bug in flexsurv for now
  }

  if(tolower(tv_type)=="baseline"){
    xbeta_m_temp <- 0
  } else if(tolower(tv_type)=="step"){
    xbeta_m_temp <- xbeta_tv_temp * as.numeric(y>tstar)
  } else if(tolower(tv_type == "logplusone")){
    xbeta_m_temp <- xbeta_tv_temp * log(y+1)
  } else if(tolower(tv_type == "piecewise")){

    # #n-vector saying which interval the ith subject falls into (when basis includes intercept)
    # cut_cats <- rowSums(basis!=0)
    # xbeta_m_temp <- xbeta_tv_temp[cbind(1:nrow(xbeta_tv_temp),cut_cats)]

    #n-vector saying which interval the ith subject falls into (when basis includes NO intercept)
    cut_cats <- rowSums(basis!=0) + 1
    #computes n length vector of the value of the time-varying x*beta associated with that interval.
    #notation is a little funky because we're indexing by a matrix with (row,column) coordinates
    #https://stackoverflow.com/questions/38271704/how-to-look-up-the-value-of-a-row-column-combination-in-a-matrix-r
    xbeta_m_temp <- cbind(0,xbeta_tv_temp)[cbind(1:nrow(xbeta_tv_temp),cut_cats)]
  }

  ll <- delta * (logh0(V_temp) - xbeta_base_temp - xbeta_m_temp) + logS0(V_temp)
  mean(-ll)
}
