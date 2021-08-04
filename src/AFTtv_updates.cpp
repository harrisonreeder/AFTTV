
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
// [[Rcpp::depends(RcppArmadillo)]]

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

