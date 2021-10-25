
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "AFTtvnew_utilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

double AFTtvnew_LN_loglik(const arma::vec &logVyL,
                            const arma::vec &logVyU,
                            const arma::vec &logVc0,
                            const arma::vec &logvyL,
                            const arma::vec &yUInf,
                            const arma::vec &yLUeq,
                            const arma::vec &c0Inf,
                            double &mu,
                            double &sigSq) {
  int n = logVyL.n_rows;
  double curr_loglik = 0;

  //in the future, I'd love to turn the yUInf etc indicators into indices so I could
  //subset vectors and use vector operations directly.

  for(int i = 0; i < n; i++){
    /* cases:
     * interval censoring (left and right diff.): compute survivor function at both, ll is log of difference.
     * (note: log of difference is also logsumexp of log survivor functions, in case it's useful)
     * right censoring (right is infinite): compute survivor function at left, ll is that
     * no censoring (left and right equal): compute log density, ll is that */
    if(yLUeq(i)==1){ //no censoring, so use log density
      curr_loglik += arma::log_normpdf(logVyL(i), mu, sqrt(sigSq))
      - logVyL(i) + logvyL(i); //note this will actually cancel out with the proposal version, but for completeness let's keep it
    } else if(yUInf(i)==1) { //right censoring, so use log survival
      curr_loglik += log1p(-arma::normcdf(logVyL(i), mu, sqrt(sigSq)));
    } else { //interval censoring, so use log of difference of survival
      curr_loglik += log(arma::normcdf(logVyU(i), mu, sqrt(sigSq))
                          - arma::normcdf(logVyL(i), mu, sqrt(sigSq)) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){ //this being 0 means the left truncation time is not log(0), so there's left truncation
      curr_loglik += -log1p(-arma::normcdf(logVc0(i), mu, sqrt(sigSq)));
    }
  }
  return curr_loglik;
}

/*BASELINE PARAMETERS (mu, sigSq)*/

void AFTtvnew_LN_update_mu(const arma::vec &logVyL,
                        const arma::vec &logVyU,
                        const arma::vec &logVc0,
                        const arma::vec &logvyL,
                        const arma::vec &yUInf,
                        const arma::vec &yLUeq,
                        const arma::vec &c0Inf,
                        double &mu,
                        double &sigSq,
                        double &mu_prop_var,
                        int &accept_mu,
                        double &curr_loglik) {
  int u;
  double loglh_prop, logR, mu_prop;
  //double logprior, logprior_prop;

  mu_prop = R::rnorm(mu, sqrt(mu_prop_var));
  loglh_prop = AFTtvnew_LN_loglik(logVyL,logVyU,logVc0,logvyL,
                         yUInf,yLUeq,c0Inf,
                         mu_prop,sigSq);

  //note, because proposal normal distribution is "symmetric"
  //aka, N(x|y,sd) = N(y|x,sd), this is just Metropolis step with no proposal distributions.
  logR = loglh_prop - curr_loglik;

  /* //if we add a prior of some kind, it should go here...
   logprior = R::dnorm(mu, mu0, sqrt(h0), 1);
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

void AFTtvnew_LN_update_sigSq(const arma::vec &logVyL,
                           const arma::vec &logVyU,
                           const arma::vec &logVc0,
                           const arma::vec &logvyL,
                           const arma::vec &yUInf,
                           const arma::vec &yLUeq,
                           const arma::vec &c0Inf,
                           double &mu,
                           double &sigSq,
                           double &sigSq_prop_var,
                           double &a_sigSq,
                           double &b_sigSq,
                           int &accept_sigSq,
                           double &curr_loglik) {
  int u;
  double loglh_prop, logR, sigSq_prop;
  double logprior, logprior_prop;

  sigSq_prop = exp( R::rnorm( log(sigSq), sqrt(sigSq_prop_var) ) );
  loglh_prop = AFTtvnew_LN_loglik(logVyL,logVyU,logVc0,logvyL,
                               yUInf,yLUeq,c0Inf,
                               mu,sigSq_prop);

  //sigSq has an inverse-gamma prior, so we define them as follows
  logprior = (-a_sigSq-1)*log(sigSq) - b_sigSq/sigSq;
  logprior_prop = (-a_sigSq-1)*log(sigSq_prop) - b_sigSq/sigSq_prop;

  //note, because proposal normal distribution is "symmetric"
  //aka, N(x|y,sd) = N(y|x,sd), this is just Metropolis step with no proposal distributions.
  //BUT, last "log(sigSq_prop) - log(sigSq)" is the "jacobian" because we sample on transformed scale
  //Because we propose from log(sigSq) to sample sigSq, so we need to add jacobian
  //https://barumpark.com/blog/2019/Jacobian-Adjustments/ has this exact example
  logR = loglh_prop - curr_loglik + logprior_prop - logprior + log(sigSq_prop) - log(sigSq);
  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    sigSq = sigSq_prop;
    accept_sigSq += 1;
    curr_loglik = loglh_prop;
  }
  return;

}

/*TIME-INVARIANT REGRESSION PARAMETERS*/

void AFTtvnew_LN_update_beta(const arma::mat &Ymat,
                          arma::vec &logVyL,
                          arma::vec &logVyU,
                          arma::vec &logVc0,
                          arma::vec &logvyL,
                          const arma::vec &yUInf,
                          const arma::vec &yLUeq,
                          const arma::vec &c0Inf,
                          const arma::mat &Xmat,
                          arma::vec &beta,
                          const arma::mat &Xvec_tv,
                          const arma::vec &beta_tv,
                          double &mu,
                          double &sigSq,
                          const arma::vec &knots,
                          double &beta_prop_var,
                          arma::vec &accept_beta,
                          double &curr_loglik){
  int p = Xmat.n_cols;
  int j, u;
  double loglh_prop, logR;
  //double logprior, logprior_prop;

  //update a single beta parameter
  arma::vec beta_prop = beta;
  j = (int) R::runif(0, p); //will need to update to sample only from non-tv X's
  beta_prop(j) = R::rnorm(beta(j), sqrt(beta_prop_var));
  //compute updated V(t) values from the proposed betaj directly
  arma::vec logVyL_prop = V_pwnew(Ymat.col(0), Xmat, beta_prop, Xvec_tv, beta_tv, knots,1);
  arma::vec logVyU_prop = V_pwnew(Ymat.col(1), Xmat, beta_prop, Xvec_tv, beta_tv, knots,1);
  arma::vec logVc0_prop = V_pwnew(Ymat.col(2), Xmat, beta_prop, Xvec_tv, beta_tv, knots,1);
  arma::vec logvyL_prop = v_pwnew(Ymat.col(0), Xmat, beta_prop, Xvec_tv, beta_tv, knots,1);

  loglh_prop = AFTtvnew_LN_loglik(logVyL_prop,logVyU_prop,logVc0_prop,logvyL_prop,
                               yUInf,yLUeq,c0Inf,
                               mu,sigSq);

  //note, because proposal normal distribution is "symmetric"
  //aka, N(x|y,sd) = N(y|x,sd), this is just Metropolis step with no proposal distributions.
  logR = loglh_prop - curr_loglik;

  /* //if we add a prior of some kind, it should go here...
   logprior
   logprior_prop
   logR += logprior_prop - logprior;
   */

  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    beta(j) = beta_prop(j);
    logVyL = logVyL_prop;
    logVyU = logVyU_prop;
    logVc0 = logVc0_prop;
    logvyL = logvyL_prop;
    accept_beta(j) += 1;
    curr_loglik = loglh_prop;
  }
  return;
}


/*TIME-VARYING REGRESSION PARAMETERS*/
void AFTtvnew_LN_update_btv(const arma::mat &Ymat,
                         arma::vec &logVyL,
                         arma::vec &logVyU,
                         arma::vec &logVc0,
                         arma::vec &logvyL,
                         const arma::vec &yUInf,
                         const arma::vec &yLUeq,
                         const arma::vec &c0Inf,
                         const arma::mat &Xmat,
                         const arma::vec &beta,
                         const arma::vec &Xvec_tv,
                         arma::vec &beta_tv,
                         double &mu,
                         double &sigSq,
                         arma::vec &knots,
                         double &btv_prop_var,
                         arma::vec &accept_btv,
                         double &curr_loglik){
  int K = beta_tv.n_rows;
  int k, u;
  double loglh_prop, logR;

  /*
  double logprior, logprior_prop;
  logprior = 0;
  logprior_prop = 0;
  */

  //choose one beta_tv element to update
  //for simplicity, we will skip the ICAR and just go with a flat prior on these for now.
  //proposal is just a normal centered at previous value for now, nothing fancy there either.
  k = (int) R::runif(0, K);
  arma::vec beta_tv_prop = beta_tv;
  beta_tv_prop(k) = R::rnorm(beta_tv(k), sqrt(btv_prop_var));

  //compute vectors with proposal V(t) values
  //later I could make this more efficient by subtracting and adding just the part that changed
  //but for now just recompute the whole thing
  arma::vec logVyL_prop = V_pwnew(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);
  arma::vec logVyU_prop = V_pwnew(Ymat.col(1), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);
  arma::vec logVc0_prop = V_pwnew(Ymat.col(2), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);

  arma::vec logvyL_prop = v_pwnew(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);

  loglh_prop = AFTtvnew_LN_loglik(logVyL_prop,logVyU_prop,logVc0_prop,logvyL_prop,
                               yUInf,yLUeq,c0Inf,
                               mu,sigSq);

  //note, because proposal normal distribution is "symmetric"
  //aka, N(x|y,sd) = N(y|x,sd), this is just Metropolis step with no proposal distributions.
  logR = loglh_prop - curr_loglik;

  //if we add a prior of some kind, it should go here...
  // here is a basic random walk prior, where the prior is centered at the previous value
  // I've just picked a default scale of 1 for the moment
  /*
   if(k>0){
   logprior = arma::log_normpdf(beta_tv(k), beta_tv(k-1), 1.0);
   logprior_prop = arma::log_normpdf(beta_tv_prop(k), beta_tv(k-1), 1.0);
   logR += logprior_prop - logprior;
   }
   */

  //Rcpp::Rcout << "beta_tvk acceptance ratio: " << exp(logR) << "\n";

  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    beta_tv = beta_tv_prop;
    logVyL = logVyL_prop;
    logVyU = logVyU_prop;
    logVc0 = logVc0_prop;
    logvyL = logvyL_prop;
    accept_btv(k) += 1;
    curr_loglik = loglh_prop;
  }
  return;
}




