
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
// [[Rcpp::depends(RcppArmadillo)]]

double AFT_LN_loglik(const arma::mat &Wmat,
                       const arma::vec &wUInf,
                       const arma::vec &wLUeq,
                       const arma::vec &c0Inf,
                       const arma::mat &Xmat,
                       arma::vec &beta,
                       double &mu,
                       double &sigSq){
  int n = Wmat.n_rows;
  double eta;
  double loglh = 0;
  arma::vec xbeta = Xmat * beta;

  //in the future, I'd love to turn the yUInf etc indicators into indices so I could
  //subset vectors and use vector operations directly.
  for(int i = 0; i < n; i++){
    eta = mu + xbeta(i);
    if(wLUeq(i)==1){ //no censoring, so use log density
      loglh      += arma::log_normpdf(Wmat(i,0), eta, sqrt(sigSq));
    } else if(wUInf(i)==1) { //right censoring, so use log survival
      loglh      += log1p(-arma::normcdf(Wmat(i,0), eta, sqrt(sigSq)));
    } else { //interval censoring, so use log of difference of survival
      loglh      += log( arma::normcdf(Wmat(i,1), eta, sqrt(sigSq))
                         - arma::normcdf(Wmat(i,0), eta, sqrt(sigSq)) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){
      loglh += -log1p(-arma::normcdf(Wmat(i,2), eta, sqrt(sigSq)));
    }
  }
  return loglh;
}




void AFT_LN_update_beta(const arma::mat &Wmat,
                 const arma::vec &wUInf,
                 const arma::vec &wLUeq,
                 const arma::vec &c0Inf,
                 const arma::mat &Xmat,
                 arma::vec &beta,
                 double &mu,
                 double &sigSq,
                 double &beta_prop_var,
                 arma::uvec &accept_beta,
                 double &curr_loglik){

  int p = Xmat.n_cols;
  int j, u;
  double loglh_prop, logR;
  //double logprior, logprior_prop;
  j = (int) R::runif(0, p);
  arma::vec beta_prop = beta;
  beta_prop(j) = R::rnorm(beta(j), sqrt(beta_prop_var));
  loglh_prop = AFT_LN_loglik(Wmat,wUInf,wLUeq,c0Inf,
                             Xmat,beta_prop,mu,sigSq);

  logR = loglh_prop - curr_loglik;

  /* //if we add a prior of some kind, it should go here...
   logprior = R::dnorm(*mu, mu0, sqrt(h0), 1);
   logprior_prop = R::dnorm(mu_prop, mu0, sqrt(h0), 1);
   logR += logprior_prop - logprior;
   */

  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    beta = beta_prop;
    accept_beta(j) += 1;
    curr_loglik = loglh_prop;
  }
  return;

}


void AFT_LN_update_mu(const arma::mat &Wmat,
               const arma::vec &wUInf,
               const arma::vec &wLUeq,
               const arma::vec &c0Inf,
               const arma::mat &Xmat,
               arma::vec &beta,
               double &mu,
               double &sigSq,
               double &mu_prop_var,
               int &accept_mu,
               double &curr_loglik){
  int u;
  //double logprior, logprior_prop;
  double mu_prop = R::rnorm(mu, sqrt(mu_prop_var));
  double loglh_prop = AFT_LN_loglik(Wmat,wUInf,wLUeq,c0Inf,
                             Xmat,beta,mu_prop,sigSq);
  double logR = loglh_prop - curr_loglik;

  /* //if we add a prior of some kind, it should go here...
   logprior = R::dnorm(*mu, mu0, sqrt(h0), 1);
   logprior_prop = R::dnorm(mu_prop, mu0, sqrt(h0), 1);
   logR += logprior_prop - logprior;
   */

  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    mu = mu_prop;
    accept_mu += 1;
    curr_loglik = loglh_prop;
  }
  return;

}

void AFT_LN_update_sigSq(const arma::mat &Wmat,
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
               int &accept_sigSq,
               double &curr_loglik){
  int u;
  double loglh_prop, logR, sigSq_prop;
  double logprior, logprior_prop;
  sigSq_prop = exp( R::rnorm( log(sigSq), sqrt(sigSq_prop_var) ) );
  loglh_prop = AFT_LN_loglik(Wmat,wUInf,wLUeq,c0Inf,
                             Xmat,beta,mu,sigSq_prop);

  //sigSq has an inverse-gamma prior, so we define them as follows
  logprior = (-a_sigSq-1)*log(sigSq) - b_sigSq/sigSq;
  logprior_prop = (-a_sigSq-1)*log(sigSq_prop) - b_sigSq/sigSq_prop;

  //Last "log(sigSq_prop) - log(sigSq)" is the "jacobian" because we sample on transformed scale

  logR = loglh_prop - curr_loglik + logprior_prop - logprior + log(sigSq_prop) - log(sigSq);
  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    sigSq = sigSq_prop;
    accept_sigSq += 1;
    curr_loglik = loglh_prop;
  }
  return;

}

