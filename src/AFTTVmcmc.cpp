
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * UPDATE FUNCTIONS
 */


void update_beta(const arma::mat &Wmat,
                 const arma::vec &wUInf,
                 const arma::vec &wLUeq,
                 const arma::vec &c0Inf,
                 const arma::mat &Xmat,
                 arma::vec &beta,
                 double &mu,
                 double &sigSq,
                 double &beta_prop_var,
                 arma::vec &accept_beta)
{

  int n = Wmat.n_rows;
  int p = Xmat.n_cols;
  int j, i, u;
  double loglh, loglh_prop, logR;
  //double logprior, logprior_prop;
  double eta, eta_prop;

  arma::vec xbeta = Xmat * beta;

  j = (int) R::runif(0, p);
  arma::vec beta_prop = beta;
  beta_prop(j) = R::rnorm(beta(j), sqrt(beta_prop_var));
  arma::vec xbeta_prop = Xmat * beta_prop;

  loglh = 0;
  loglh_prop = 0;

  for(i = 0; i < n; i++){
    eta = mu + xbeta(i);
    eta_prop = mu + xbeta_prop(i);
    if(wLUeq(i)==1){ //no censoring, so use log density
      loglh      += arma::log_normpdf(Wmat(i,0), eta, sqrt(sigSq));
      loglh_prop += arma::log_normpdf(Wmat(i,0), eta_prop, sqrt(sigSq));
      // loglh      += R::dnorm(Wmat(i,0), eta, sqrt(sigSq), 1);
      // loglh_prop += R::dnorm(Wmat(i,0), eta_prop, sqrt(sigSq), 1);
    } else if(wUInf(i)==1) { //right censoring, so use log survival
      loglh      += log1p(-arma::normcdf(Wmat(i,0), eta, sqrt(sigSq)));
      loglh_prop += log1p(-arma::normcdf(Wmat(i,0), eta_prop, sqrt(sigSq)));
      // loglh      += R::pnorm(Wmat(i,0), eta, sqrt(sigSq), 0, 1);
      // loglh_prop += R::pnorm(Wmat(i,0), eta_prop, sqrt(sigSq), 0, 1);
    } else { //interval censoring, so use log of difference of survival
      loglh      += log( arma::normcdf(Wmat(i,1), eta, sqrt(sigSq))
                           - arma::normcdf(Wmat(i,0), eta, sqrt(sigSq)) );
      loglh_prop += log( arma::normcdf(Wmat(i,1), eta_prop, sqrt(sigSq))
                           - arma::normcdf(Wmat(i,0), eta_prop, sqrt(sigSq)) );
      // loglh      += log( R::pnorm(Wmat(i,0), eta, sqrt(sigSq), 0, 0)
      //                      - R::pnorm(Wmat(i,1), eta, sqrt(sigSq), 0, 0) );
      // loglh_prop += log( R::pnorm(Wmat(i,0), eta_prop, sqrt(sigSq), 0, 0)
      //                      - R::pnorm(Wmat(i,1), eta_prop, sqrt(sigSq), 0, 0) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){
      loglh += -log1p(-arma::normcdf(Wmat(i,2), eta, sqrt(sigSq)));
      loglh_prop += -log1p(-arma::normcdf(Wmat(i,2), eta_prop, sqrt(sigSq)));
      // loglh += -R::pnorm(Wmat(i,2), eta, sqrt(sigSq), 0, 1);
      // loglh_prop += -R::pnorm(Wmat(i,2), eta_prop, sqrt(sigSq), 0, 1);
    }
  }
  logR = loglh_prop - loglh;

  /* //if we add a prior of some kind, it should go here...
   logprior = R::dnorm(*mu, mu0, sqrt(h0), 1);
   logprior_prop = R::dnorm(mu_prop, mu0, sqrt(h0), 1);
   logR += logprior_prop - logprior;
   */

  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    beta = beta_prop;
    accept_beta(j) += 1;
  }
  return;

}


void update_mu(const arma::mat &Wmat,
               const arma::vec &wUInf,
               const arma::vec &wLUeq,
               const arma::vec &c0Inf,
               const arma::mat &Xmat,
               arma::vec &beta,
               double &mu,
               double &sigSq,
               double &mu_prop_var,
               int &accept_mu)
{
  int n = Wmat.n_rows;
  int i, u;
  double loglh, loglh_prop, logR, mu_prop;
  //double logprior, logprior_prop;
  double eta, eta_prop;
  arma::vec xbeta = Xmat * beta;

  loglh = 0;
  loglh_prop = 0;
  mu_prop = R::rnorm(mu, sqrt(mu_prop_var));

  for(i = 0; i < n; i++){
    eta = mu + xbeta(i);
    eta_prop = mu_prop + xbeta(i);
    /* cases:
     * interval censoring (left and right diff.): compute survivor function at both, ll is log of difference.
     * (note: log of difference is also logsumexp of log survivor functions, in case it's useful)
     * right censoring (right is infinite): compute survivor function at left, ll is that
     * no censoring (left and right equal): compute log density, ll is that */
    if(wLUeq(i)==1){ //no censoring, so use log density
      loglh      += arma::log_normpdf(Wmat(i,0), eta, sqrt(sigSq));
      loglh_prop += arma::log_normpdf(Wmat(i,0), eta_prop, sqrt(sigSq));
      // loglh      += R::dnorm(Wmat(i,0), eta, sqrt(sigSq), 1);
      // loglh_prop += R::dnorm(Wmat(i,0), eta_prop, sqrt(sigSq), 1);
    } else if(wUInf(i)==1) { //right censoring, so use log survival
      loglh      += log1p(-arma::normcdf(Wmat(i,0), eta, sqrt(sigSq)));
      loglh_prop += log1p(-arma::normcdf(Wmat(i,0), eta_prop, sqrt(sigSq)));
      // loglh      += R::pnorm(Wmat(i,0), eta, sqrt(sigSq), 0, 1);
      // loglh_prop += R::pnorm(Wmat(i,0), eta_prop, sqrt(sigSq), 0, 1);
    } else { //interval censoring, so use log of difference of survival
      loglh      += log( arma::normcdf(Wmat(i,1), eta, sqrt(sigSq))
                           - arma::normcdf(Wmat(i,0), eta, sqrt(sigSq)) );
      loglh_prop += log( arma::normcdf(Wmat(i,1), eta_prop, sqrt(sigSq))
                           - arma::normcdf(Wmat(i,0), eta_prop, sqrt(sigSq)) );
      // loglh      += log( R::pnorm(Wmat(i,0), eta, sqrt(sigSq), 0, 0)
      //                      - R::pnorm(Wmat(i,1), eta, sqrt(sigSq), 0, 0) );
      // loglh_prop += log( R::pnorm(Wmat(i,0), eta_prop, sqrt(sigSq), 0, 0)
      //                      - R::pnorm(Wmat(i,1), eta_prop, sqrt(sigSq), 0, 0) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){
      loglh += -log1p(-arma::normcdf(Wmat(i,2), eta, sqrt(sigSq)));
      loglh_prop += -log1p(-arma::normcdf(Wmat(i,2), eta_prop, sqrt(sigSq)));
      // loglh += -R::pnorm(Wmat(i,2), eta, sqrt(sigSq), 0, 1);
      // loglh_prop += -R::pnorm(Wmat(i,2), eta_prop, sqrt(sigSq), 0, 1);
    }
  }
  logR = loglh_prop - loglh;

  /* //if we add a prior of some kind, it should go here...
   logprior = R::dnorm(*mu, mu0, sqrt(h0), 1);
   logprior_prop = R::dnorm(mu_prop, mu0, sqrt(h0), 1);
   logR += logprior_prop - logprior;
   */

  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    mu = mu_prop;
    accept_mu += 1;
  }
  return;

}

void update_sigSq(const arma::mat &Wmat,
               const arma::vec &wUInf,
               const arma::vec &wLUeq,
               const arma::vec &c0Inf,
               const arma::mat &Xmat,
               arma::vec &beta,
               double &mu,
               double &sigSq,
               double &sigSq_prop_var,
               double &a_sigSq,
               double &b_sigSq,
               int &accept_sigSq)
{
  int n = Wmat.n_rows;
  int i, u;
  double loglh, loglh_prop, logR, sigSq_prop;
  double logprior, logprior_prop;
  double eta;
  arma::vec xbeta = Xmat * beta;

  loglh = 0;
  loglh_prop = 0;
  sigSq_prop = exp( R::rnorm( log(sigSq), sqrt(sigSq_prop_var) ) );

  for(i = 0; i < n; i++){
    eta = mu + xbeta(i);
    if(wLUeq(i)==1){ //no censoring, so use log density
      loglh      += arma::log_normpdf(Wmat(i,0), eta, sqrt(sigSq));
      loglh_prop += arma::log_normpdf(Wmat(i,0), eta, sqrt(sigSq_prop));
      // loglh      += R::dnorm(Wmat(i,0), eta, sqrt(sigSq), 1);
      // loglh_prop += R::dnorm(Wmat(i,0), eta, sqrt(sigSq_prop), 1);
    } else if(wUInf(i)==1) { //right censoring, so use log survival
      loglh      += log1p(-arma::normcdf(Wmat(i,0), eta, sqrt(sigSq)));
      loglh_prop += log1p(-arma::normcdf(Wmat(i,0), eta, sqrt(sigSq_prop)));
      // loglh      += R::pnorm(Wmat(i,0), eta, sqrt(sigSq), 0, 1);
      // loglh_prop += R::pnorm(Wmat(i,0), eta, sqrt(sigSq_prop), 0, 1);
    } else { //interval censoring, so use log of difference of survival
      loglh      += log( arma::normcdf(Wmat(i,1), eta, sqrt(sigSq))
                           - arma::normcdf(Wmat(i,0), eta, sqrt(sigSq)) );
      loglh_prop += log( arma::normcdf(Wmat(i,1), eta, sqrt(sigSq_prop))
                           - arma::normcdf(Wmat(i,0), eta, sqrt(sigSq_prop)) );
      // loglh      += log( R::pnorm(Wmat(i,0), eta, sqrt(sigSq), 0, 0)
      //                    - R::pnorm(Wmat(i,1), eta, sqrt(sigSq), 0, 0) );
      // loglh_prop += log( R::pnorm(Wmat(i,0), eta, sqrt(sigSq_prop), 0, 0)
      //                    - R::pnorm(Wmat(i,1), eta, sqrt(sigSq_prop), 0, 0) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){
      loglh += -log1p(-arma::normcdf(Wmat(i,2), eta, sqrt(sigSq)));
      loglh_prop += -log1p(-arma::normcdf(Wmat(i,2), eta, sqrt(sigSq_prop)));
      // loglh += -R::pnorm(Wmat(i,2), eta, sqrt(sigSq), 0, 1);
      // loglh_prop += -R::pnorm(Wmat(i,2), eta, sqrt(sigSq_prop), 0, 1);
    }
  }

  //sigSq has an inverse-gamma prior, so we define them as follows
  logprior = (-a_sigSq-1)*log(sigSq) - b_sigSq/sigSq;
  logprior_prop = (-a_sigSq-1)*log(sigSq_prop) - b_sigSq/sigSq_prop;

  //Last "log(sigSq_prop) - log(sigSq)" is the "jacobian" because we sample on transformed scale

  logR = loglh_prop - loglh + logprior_prop - logprior + log(sigSq_prop) - log(sigSq);
  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    sigSq = sigSq_prop;
    accept_sigSq += 1;
  }
  return;

}



/*
 * SAMPLER FUNCTION
 */

// [[Rcpp::export]]
Rcpp::List BAFTtvLTmcmc(const arma::mat& Wmat,
                        const arma::vec& wUInf,
                        const arma::vec& wLUeq,
                        const arma::vec& c0Inf,
                        const arma::mat& Xmat,
                        const arma::vec& hyperP,
                        const arma::vec& mcmcP,
                        const arma::vec& startValues,
                        int n_burnin,
                        int n_sample,
                        int thin){
  //timekeeping objects
  std::time_t newt;

  //set constants
  int n = Xmat.n_rows;
  int p = Xmat.n_cols;
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'
  int move; //index for which parameter to update
  int M; //counter for MCMC sampler
  int StoreInx; //index for where to store a sample, post-thinning

  //set hyperparameters
  double a_sigSq = hyperP(0);
  double b_sigSq = hyperP(1);

  //set MCMC tuning parameters
  double beta_prop_var = mcmcP(0);
  double mu_prop_var = mcmcP(1);
  double sigSq_prop_var = mcmcP(2);

  /*
  double pBeta = 1.0/3.0; //probability of updating beta (random scan)
  double pMu = 1.0/3.0; //probability of updating mu (random scan)
  double pSigSq = 1.0/3.0;  //probability of updating sigSq (random scan)

  arma::vec moveProb = arma::vec(3);
  moveProb(0) = pBeta;
  moveProb(1) = pMu;
  moveProb(2) = pSigSq;
  */


  //initialize starting values
  double mu = startValues(0);
  double sigSq = startValues(1);
  arma::vec beta = startValues(arma::span(2,2+p-1));

  //create storage objects
  arma::mat sample_beta = arma::mat(p,n_store); //store betas column by column for "speed" ?
  arma::vec sample_mu = arma::vec(n_store);
  arma::vec sample_sigSq = arma::vec(n_store);
  arma::vec accept_beta = arma::vec(p);
  int accept_mu = 0;
  int accept_sigSq = 0;

  //eventually, compute V(t) but for now stick with W

  for(M = 0; M < n_iter; M++){

    //if we've changed beta, recompute V(t) (or in time-invariant case just eta), otherwise no need to

    /*
    update_beta(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, beta_prop_var, accept_beta);
    update_mu(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, mu_prop_var, accept_mu);
    update_sigSq(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, sigSq_prop_var, a_sigSq, b_sigSq, accept_sigSq);
    */

    move = (int) R::runif(0, 3); //for now, equal probability of each of the 3 moves.

    if(move == 0){
      update_beta(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, beta_prop_var, accept_beta);
    }
    if(move == 1){
      update_mu(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, mu_prop_var, accept_mu);
    }
    if(move == 2){
      update_sigSq(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, sigSq_prop_var, a_sigSq, b_sigSq, accept_sigSq);
    }

    /* Storing posterior samples */
    if( ( (M+1) % thin ) == 0 && (M+1) > n_burnin)
    {
      StoreInx = (M+1 - n_burnin)/thin;

      sample_beta.col(StoreInx - 1) = beta;
      sample_mu(StoreInx - 1) = mu;
      sample_sigSq(StoreInx - 1) = sigSq;
    }


    if( ( (M+1) % 10000 ) == 0){
      newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      Rcpp::Rcout << "iteration: " << M+1 << ": " << ctime(&newt) << "\n";
      Rcpp::checkUserInterrupt(); //cheks if the user hit the "stop" icon to cancel running sampler.

      //this is this like alternative
      // Rcpp::Rcout << "iteration: " << M+1 << ": " << ctime_s(tmBuff, sizeof(tmBuff), &newt) << "\n";
    }

  }

  return Rcpp::List::create(
    Rcpp::Named("samples") = Rcpp::List::create(
      Rcpp::Named("beta") = sample_beta.t(),
      Rcpp::Named("mu") = sample_mu,
      Rcpp::Named("sigSq") = sample_sigSq),
    Rcpp::Named("accept") = Rcpp::List::create(
      Rcpp::Named("beta") = accept_beta,
      Rcpp::Named("mu") = accept_mu,
      Rcpp::Named("sigSq") = accept_sigSq)
  );

}

