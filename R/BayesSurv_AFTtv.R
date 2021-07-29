
#' Bayesian Accelerated Failure Time Model with Percentile-varying Effects
#'
#' @param Y a matrix with three columns: left end of interval, right end of interval, left truncation time
#' @param Xmat
#' @param hyperParams
#' @param startValues
#' @param mcmcParams
#'
#' @return
#' @export
BayesSurv_AFTtv <- function(Y,
                     Xmat,
                     hyperParams,
                     startValues_vec,
                     mcmcParams)
{
  # browser()
  #this version only fits one chain, but we can get multiple chains
  #by externally running twice?
  ###
  n	<- dim(Y)[1]
  p	<- ncol(Xmat)

  ###
  # hyperP  <- as.vector(c(hyperParams$LN$a.sigSq, hyperParams$LN$b.sigSq#,
  #                        # hyperParams$LN$mu0,hyperParams$LN$h0
  #                        ))
  hyperP <- hyperParams$LN$LN.ab
  #NOTE I'M RETROFITTING THE WORD 'ZETA' HERE BUT WILL CHANGE
  mcmcP   <- as.vector(c(mcmcParams$tuning$beta.prop.var,
                         mcmcParams$tuning$mu.prop.var,
                         mcmcParams$tuning$zeta.prop.var)) #should be sigSq

  numReps     <- mcmcParams$run$numReps
  thin        <- mcmcParams$run$thin
  burninPerc  <- mcmcParams$run$burninPerc
  n_burnin <- numReps * burninPerc
  n_sample <- numReps - n_burnin
  n_store      <- numReps/thin * (1 - burninPerc)

  stopifnot(n_burnin %% 1 == 0)
  stopifnot(n_sample > 0)
  if(n_store %% 1 != 0){
    stop("numReps * burninPerc  must be divisible by thin")
  }

  W <- Y
  W[,1] <- log(Y[,1])
  W[,2] <- log(Y[,2])
  W[,3] <- log(Y[,3])

  for(i in 1:n) if(W[i,1] == -Inf)
  {
    W[i,1] <- -9.9e10
  }

  wUInf <- rep(0, n)
  for(i in 1:n) if(W[i,2] == Inf)
  {
    W[i,2] <- 9.9e10
    wUInf[i] <- 1
  }

  c0Inf <- rep(0, n)
  for(i in 1:n) if(W[i,3] == -Inf)
  {
    W[i,3] <- -9.9e10
    c0Inf[i] <- 1
  }

  mcmcRet     <- BAFTtvLTmcmc(
                    Wmat        = W,
                    wUInf			  = wUInf,
                    c0Inf			  = c0Inf,
                    Xmat        = Xmat,
                    hyperP      = hyperP,
                    mcmcP       = mcmcP,
                    startValues = startValues_vec,
                    n_burnin		= numReps,
                    n_sample		= n_sample,
                    thin        = thin)
  return(mcmcRet)

  # w.p <- matrix(as.vector(mcmcRet$samples_w), nrow=nStore, byrow=T)
  # if(p >0){
  #   beta.p <- matrix(as.vector(mcmcRet$samples_beta), nrow=nStore, byrow=T)
  # } else{
  #   beta.p <- NULL
  # }
  #
  # mu.p <- matrix(as.vector(mcmcRet$samples_mu), nrow=nStore, byrow=T)
  # sigSq.p <- matrix(as.vector(mcmcRet$samples_sigSq), nrow=nStore, byrow=T)
  # lambdaSq.p <- matrix(as.vector(mcmcRet$samples_lambdaSq), nrow=nStore, byrow=T)
  #
  # if(p >0)
  # {
  #   accept.beta	 <- as.vector(mcmcRet$samples_misc[1:p])
  # }else
  # {
  #   accept.beta <- NULL
  # }
  #
  # accept.mu	 <- as.vector(mcmcRet$samples_misc[p+1])
  # accept.sigSq	 <- as.vector(mcmcRet$samples_misc[p+2])
  #
  # ret <- list(beta.p = beta.p, mu.p=mu.p, sigSq.p = sigSq.p,
  #             accept.beta = accept.beta, accept.mu = accept.mu, accept.sigSq = accept.sigSq)
  #
  # class(ret) <- "BAFTtv"
  # return(ret)
}
