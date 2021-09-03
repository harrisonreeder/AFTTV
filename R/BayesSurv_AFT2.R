
#' Bayesian Accelerated Failure Time Model
#'
#' @param Y an n by 3 matrix with columns: left end of interval, right end of interval, left truncation time
#' @param Xmat an n by p matrix of covariates
#' @param Xvec_tv an n-length vector with the covariate having a time-varying effect
#' @param prior_list a list of priors being used
#' @param tuning_vec a vec of tuning parameters
#' @param hyper_list a list of hyperparameters
#' @param knots a vector of knots (doesn't start with 0, does end with maximum time)
#' @param numReps numeric total number of iterations
#' @param thin numeric number of iterations to take before each sample
#' @param burninPerc numeric proportion of numReps to toss as burn-in
#' @param n_chains numeric number of chains to initialize, if start_list is null
#' @param start_list a list of lists of start values. If null, n_chains sets number of chains
#'
#' @return a list with outputs
#' @export
BayesSurv_AFT2 <- function(Y,
                           Xmat,
                           # prior_list,
                           tuning_vec,
                           # tuning_list,
                           hyper_list,
                           numReps,
                           thin,
                           burninPerc,
                           n_chains=NULL,
                           start_list=NULL){
  # browser()
  #set index values
  n	<- dim(Y)[1]
  p	<- ncol(Xmat) #TODO: figure out what to do if no covariates

  #set mcmc parameters
  n_burnin <- numReps * burninPerc
  n_sample <- numReps - n_burnin
  n_store      <- numReps/thin * (1 - burninPerc)
  stopifnot(n_burnin %% 1 == 0)
  stopifnot(n_sample > 0)
  if(n_store %% 1 != 0){ stop("numReps * burninPerc  must be divisible by thin")}


  # ####SET PRIOR SPECIFICATIONS####
  # #convert prior_list into list of numbers that can be passed to cpp
  # #0 is flat prior, 1 is inv-gamma prior, 2 is mvn-icar prior
  # prior_list_num <- prior_list
  # prior_list_num <- rapply(prior_list_num,function(x) ifelse(tolower(x) %in% c("flat"),0,x), how = "replace")
  # prior_list_num <- rapply(prior_list_num,function(x) ifelse(tolower(x) %in% c("invgamma","ig","inv-gamma"),1,x), how = "replace")
  # #add automatic filling in of defaults
  # stopifnot(all(names(prior_list_num) %in% c("mu","sigSq","beta"))) #these all need SOMETHING
  # prior_vec_num <- unlist(prior_list_num[c("mu","sigSq","beta")])

  # #instead of using the prior_list, we could also build from a prior_vec
  # prior_vec_num <- sapply(tolower(prior_vec), switch,
  #                         "flat"=0,"invgamma"=1,"inv-gamma"=1,
  #                         "mvn-icar"=2,"mvnicar"=2,"icar"=2)


  ####ASSIGN START VALUES####
  if(is.null(start_list)){
    if(is.null(n_chains)){stop("either supply non-null start_list, or specify n_chains.")}
    #create an empty list with as many empty lists as there are chains
    start_list <- list()
    for(i in 1:n_chains){start_list[[i]] <- list()}
  }
  #loop through and fill in any start values not already initialized
  for(i in 1:n_chains){
    #TODO: what if no betas are specified?

    #these randomization choices match kyu ha
    if(is.null(start_list[[i]][["mu"]])){ start_list[[i]][["mu"]] <- runif(1, -0.1, 0.1) }
    if(is.null(start_list[[i]][["sigSq"]])){ start_list[[i]][["sigSq"]] <- runif(1, 0.5, 1.5) }
    if(is.null(start_list[[i]][["beta"]])){ start_list[[i]][["beta"]] <- runif(p, -0.1, 0.1) }
    #ensure everything is in the right order
    start_list[[i]][c("mu","sigSq","beta")]
  }
  start_mat <- sapply(start_list,unlist)

  ####CHECK HYPERPARAMETERS####
  # if(prior_list_num$sigSq==1){ stopifnot(length(hyper_list$sigSq)==2)}
  hyper_vec <- unlist(hyper_list["sigSq"])

  ####PREP DATA####
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

  wLUeq <- as.numeric(W[,1] == W[,2])

  c0Inf <- rep(0, n)
  for(i in 1:n) if(W[i,3] == -Inf)
  {
    W[i,3] <- -9.9e10
    c0Inf[i] <- 1
  }

  # #total number of parameters sampled
  p_tot <- length(unlist(start_list[[1]]))
  # #TODO: parallelize this loop (I know it can be done!)
  out_list <- list()
  # #generate an array to store the resulting samples
  out_list[["samples"]] <- array(dim = c(n_store, n_chains, p_tot))
  out_list[["accept"]] <- list()
  for(i in 1:n_chains){
    print(paste0("Chain: ", i))
    #mcmcRet[[i]]     <- AFTtv_LN_mcmc( #for now, set this aside so I can test my rj version
    mcmcRet     <- AFT_LN_mcmc(
      Wmat        = W,
      wUInf			  = wUInf,
      wLUeq			  = wLUeq,
      c0Inf			  = c0Inf,
      Xmat        = Xmat,
      hyper_vec      = hyper_vec,
      tuning_vec       = tuning_vec,
      start_vec = start_mat[,i],
      n_burnin		= n_burnin,
      n_sample		= n_sample,
      thin        = thin)
    out_list[["samples"]][,i,] <- do.call(what = cbind,
                                          args = mcmcRet$samples)
    out_list[["accept"]][[paste0("chain",i)]] <- mcmcRet$accept
  }

  dimnames(out_list[["samples"]]) <- list(as.character(1:n_store),
                                          paste0("chain:",1:n_chains),
                                          c("mu","sigSq",paste0("beta",1:p)))

  out_list[["setup"]]	<-
    list(hyper_list = hyper_list,
         start_list = start_list,
         prior_list = prior_list,
         mcmc_para = c(n_burnin=n_burnin,n_sample=n_sample,thin=thin))

  return(out_list)


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

