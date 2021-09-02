
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "AFTtv_utilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

double AFTtv_LN_loglik(const arma::vec &logVyL,
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

void AFTtv_LN_update_mu(const arma::vec &logVyL,
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
  int n = logVyL.n_rows;
  int u;
  double loglh_prop, logR, mu_prop;
  //double logprior, logprior_prop;

  mu_prop = R::rnorm(mu, sqrt(mu_prop_var));
  loglh_prop = AFTtv_LN_loglik(logVyL,logVyU,logVc0,logvyL,
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

void AFTtv_LN_update_sigSq(const arma::vec &logVyL,
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
  int n = logVyL.n_rows;
  int u;
  double loglh_prop, logR, sigSq_prop;
  double logprior, logprior_prop;

  sigSq_prop = exp( R::rnorm( log(sigSq), sqrt(sigSq_prop_var) ) );
  loglh_prop = AFTtv_LN_loglik(logVyL,logVyU,logVc0,logvyL,
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

void AFTtv_LN_update_beta(arma::vec &logVyL,
                          arma::vec &logVyU,
                          arma::vec &logVc0,
                          arma::vec &logvyL,
                          const arma::vec &yUInf,
                          const arma::vec &yLUeq,
                          const arma::vec &c0Inf,
                          const arma::mat &Xmat,
                          arma::vec &beta,
                          double &mu,
                          double &sigSq,
                          double &beta_prop_var,
                          arma::vec &accept_beta,
                          double &curr_loglik){
  int n = logVyL.n_rows;
  int p = Xmat.n_cols;
  int j, u;
  double loglh_prop, logR;
  //double logprior, logprior_prop;

  //update a single beta parameter
  j = (int) R::runif(0, p); //will need to update to sample only from non-tv X's
  double betaj_prop = R::rnorm(beta(j), sqrt(beta_prop_var));
  //compute updated V(t) values from the proposed betaj directly
  arma::vec xbetaj_propdiff = Xmat.col(j) * (beta(j) - betaj_prop);
  arma::vec logVyL_prop = logVyL + xbetaj_propdiff;
  arma::vec logVyU_prop = logVyU + xbetaj_propdiff;
  arma::vec logVc0_prop = logVc0 + xbetaj_propdiff;
  arma::vec logvyL_prop = logvyL + xbetaj_propdiff;

  loglh_prop = AFTtv_LN_loglik(logVyL_prop,logVyU_prop,logVc0_prop,logvyL_prop,
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
    beta(j) = betaj_prop;
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

void AFTtv_LN_update_hyper_icar(const arma::vec &beta_tv,
                                double &meanbtv,
                                double &varbtv,
                                double &a_btv,
                                double &b_btv,
                                const arma::mat &cholinvSigma_btv){

  int K = beta_tv.n_rows; //update this once beta can change dimension
  double mean_num, var_num, denom;
  double meanbtv_mean, meanbtv_var;
  double a_part, b_part;

  //prepare some quadratic forms used to define the conjugate posterior distros for hyperparams
  arma::vec Rones = arma::ones(K);
  inplace_tri_mat_mult(cholinvSigma_btv(arma::span(0,K-1),arma::span(0,K-1)), Rones);
  //sets Rones to 1.t() * R.t() * R * 1 = 1.t() * invSigma * 1
  denom = arma::dot(Rones,Rones);
  arma::vec Rbeta_tv = beta_tv;
  inplace_tri_mat_mult(cholinvSigma_btv(arma::span(0,K-1),arma::span(0,K-1)), Rbeta_tv);
  //sets Rbeta_tv to 1.t() * R.t() * R * beta_tv = 1.t() * invSigma * beta_tv
  mean_num = arma::dot(Rones,Rbeta_tv);

  //Gibbs sample new mean hyperparameter
  meanbtv_mean = mean_num/denom;
  meanbtv_var = varbtv/denom;
  //Gibbs draw from conjugate posterior
  meanbtv = R::rnorm(meanbtv_mean, sqrt(meanbtv_var));


  //computes (meanbtv*1 - beta_tv).t() * invSigma * (meanbtv*1 - beta_tv)
  a_part = a_btv + K/2.0;
  b_part = b_btv + arma::dot( (meanbtv * Rones - Rbeta_tv),
                              (meanbtv * Rones - Rbeta_tv) )/2;
  //Kyu ha's notes have the "b" parameter in terms of rate, but we invert to put in terms of scale.
  //we draw with gibbs, then invert the drawn parameter because the conjugacy is with inverse of the var param.
  varbtv = 1/R::rgamma(a_part,1/b_part);

  return;

}


void AFTtv_LN_update_btv_icar(const arma::mat &Ymat,
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
                         double &meanbtv,
                         double &varbtv,
                         const arma::mat &cholinvSigma_btv,
                         arma::vec &accept_btv,
                         double &curr_loglik){

  int n = logVyL.n_rows;
  int K = beta_tv.n_rows; //update this once beta_tv has changing dimension
  int k, u;
  double loglh_prop, logprior, logprior_prop, logR;
//  double nu_btv, nu_btv_prop;

  logprior = 0;
  logprior_prop = 0;

  //choose one beta_tv element to update
  //proposal is just a normal centered at previous value for now, nothing fancy with newton steps;
  k = (int) R::runif(0, K);
  arma::vec beta_tv_prop = beta_tv;
  beta_tv_prop(k) = R::rnorm(beta_tv(k), sqrt(btv_prop_var));

  //compute vectors with proposal V(t) values
  //later I could make this more efficient by subtracting and adding just the part that changed
  //but for now just recompute the whole thing
  arma::vec logVyL_prop = Vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);
  arma::vec logVyU_prop = Vx_pw(Ymat.col(1), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);
  arma::vec logVc0_prop = Vx_pw(Ymat.col(2), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);

  arma::vec logvyL_prop = vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);

  loglh_prop = AFTtv_LN_loglik(logVyL_prop,logVyU_prop,logVc0_prop,logvyL_prop,
                               yUInf,yLUeq,c0Inf,
                               mu,sigSq);

  //compute prior information for MVN-ICAR
  //If we don't do newton-raphson then it's pretty straightforward actually!
  if(K>1){ //if there are multiple intervals
    logprior = dmvnrm_arma(beta_tv, meanbtv, varbtv,
                           cholinvSigma_btv(arma::span(0,K-1),arma::span(0,K-1)),1);
    logprior_prop = dmvnrm_arma(beta_tv_prop, meanbtv, varbtv,
                                cholinvSigma_btv(arma::span(0,K-1),arma::span(0,K-1)),1);
  } else{ //if there's just one interval
    //here, we just feed in the first (and only) value of beta_tv,
    //(note, in the univariate case knots(0) IS the length of the interval, and therefore the value of Q(0) = Sigma(0,0)
    logprior = arma::log_normpdf(beta_tv(0), meanbtv, sqrt(varbtv*knots(0)));
    logprior_prop = arma::log_normpdf(beta_tv_prop(0), meanbtv, sqrt(varbtv*knots(0)));
  }

  //note, because for now proposal normal distribution is "symmetric"
  //aka, N(x|y,sd) = N(y|x,sd), this is just Metropolis step with no proposal distributions.
  //no need for any other jacobian corrections because we are not proposing on a different scale
  logR = loglh_prop - curr_loglik + logprior_prop - logprior;

  /* //Don't need to compute nu because it's only used for newton steps
  if(K > 1) {
    if(k == 0){
      nu_btv = nu_btv - IminusW(0,1) * (beta_tv(1) - mu_btv);
    } else if(k == K-1){
      nu_btv = nu_btv - IminusW(K-1,K-2) * (beta_tv(K-2) - mu_btv);
    } else{
      nu_btv = nu_btv - IminusW(k,k-1) * (beta_tv(k-1) - mu_btv)
                      - IminusW(k,k+1) * (beta_tv(k+1) - mu_btv);
    }
  } else{
    nu_btv = mu_btv;
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


void AFTtv_LN_update_btv(const arma::mat &Ymat,
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

  int n = logVyL.n_rows;
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
  arma::vec logVyL_prop = Vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);
  arma::vec logVyU_prop = Vx_pw(Ymat.col(1), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);
  arma::vec logVc0_prop = Vx_pw(Ymat.col(2), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);

  arma::vec logvyL_prop = vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);

  loglh_prop = AFTtv_LN_loglik(logVyL_prop,logVyU_prop,logVc0_prop,logvyL_prop,
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




