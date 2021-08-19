
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "AFTtv_LN_updates.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List AFTtv_LN_mcmc(const arma::mat &Ymat,
                        const arma::vec &yUInf,
                        const arma::vec &yLUeq,
                        const arma::vec &c0Inf,
                        const arma::mat &Xmat,
                        const arma::vec &Xvec_tv,
                        const arma::vec &hyperP,
                        double beta_prop_var,
                        double btv_prop_var,
                        double mu_prop_var,
                        double sigSq_prop_var,
                        const arma::vec &knots_init,
                        const arma::vec &startValues,
                        int n_burnin,
                        int n_sample,
                        int thin){
  //timekeeping objects
  std::time_t newt;

  //set constants
  int p = Xmat.n_cols;
  //one parameter for every element of the "knots" vector, because it excludes 0 but includes infty.
  int p_tv = knots_init.n_rows;
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'
  int move; //index for which parameter to update
  int M; //counter for MCMC sampler
  int StoreInx; //index for where to store a sample, post-thinning

  //set hyperparameters
  double a_sigSq = hyperP(0);
  double b_sigSq = hyperP(1);
//  double a_btv = hyperP(2); //prior marginal mean for time-varying beta coefficients
//  double b_btv = hyperP(3); //prior variance for time-varying beta coefficients
//  double c_btv = hyperP(4); //prior dependence parameter for time-varying beta coefficients

  /* manually set random scan probabilities
  double pBetaTV = 1.0/4.0; //probability of updating beta (random scan)
  double pBeta = 1.0/4.0; //probability of updating beta (random scan)
  double pMu = 1.0/4.0; //probability of updating mu (random scan)
  double pSigSq = 1.0/4.0;  //probability of updating sigSq (random scan)

  arma::vec moveProb = arma::vec(4);
  moveProb(0) = pBeta;
  moveProb(1) = pBetaTV;
  moveProb(2) = pMu;
  moveProb(3) = pSigSq;
  */

  //initialize starting values
  double mu = startValues(0);
  double sigSq = startValues(1);
  arma::vec beta = startValues(arma::span(2,2+p-1));
  arma::vec beta_tv = startValues(arma::span(2+p,2+p-1+p_tv));
  arma::vec knots = knots_init;

  //create storage objects for results
  arma::mat sample_beta = arma::mat(p,n_store); //store betas column by column for "speed" ?
  arma::mat sample_btv = arma::mat(p_tv,n_store); //store betas column by column for "speed" ?
  arma::vec sample_mu = arma::vec(n_store);
  arma::vec sample_sigSq = arma::vec(n_store);
  arma::vec accept_beta = arma::vec(p);
  arma::vec accept_btv = arma::vec(p_tv);
  int accept_mu = 0;
  int accept_sigSq = 0;

  //initialize V(t) values for yL, yU and cL
  arma::vec logVyL = Vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  arma::vec logVyU = Vx_pw(Ymat.col(1), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  arma::vec logVc0 = Vx_pw(Ymat.col(2), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  //initialize v(t) values for yL (only used when observation is exactly observed)
  arma::vec logvyL = vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv, knots, 1);

//  Rcpp::Rcout << "initial logVyL: " << logVyL(arma::span(0,10)) << "\n";

  for(M = 0; M < n_iter; M++){
//    Rcpp::Rcout << "iteration: " << M << "\n";

    AFTtv_LN_update_mu(logVyL, logVyU, logVc0, logvyL,
                       yUInf, yLUeq, c0Inf,
                       mu, sigSq, mu_prop_var, accept_mu);
//    Rcpp::Rcout << "updated mu: " << mu << "\n";
    AFTtv_LN_update_sigSq(logVyL, logVyU, logVc0, logvyL,
                          yUInf, yLUeq, c0Inf,
                          mu, sigSq, sigSq_prop_var,
                          a_sigSq, b_sigSq, accept_sigSq);
//    Rcpp::Rcout << "updated sigSq: " << sigSq << "\n";
    AFTtv_LN_update_beta(logVyL, logVyU, logVc0, logvyL,
                         yUInf, yLUeq, c0Inf,
                         Xmat, beta, mu, sigSq, beta_prop_var, accept_beta);
//    Rcpp::Rcout << "updated beta: " << beta.t() << "\n";
    AFTtv_LN_update_btv(Ymat, logVyL, logVyU, logVc0, logvyL,
                         yUInf, yLUeq, c0Inf,
                         Xmat, beta, Xvec_tv, beta_tv,
                         mu, sigSq, knots, btv_prop_var, accept_btv);
//    Rcpp::Rcout << "updated beta_tv: " << beta_tv.t() << "\n";

    // Storing posterior samples
    if( ( (M+1) % thin ) == 0 && (M+1) > n_burnin)
    {
      StoreInx = (M+1 - n_burnin)/thin;

      sample_beta.col(StoreInx - 1) = beta;
      sample_btv.col(StoreInx - 1) = beta_tv;
      sample_mu(StoreInx - 1) = mu;
      sample_sigSq(StoreInx - 1) = sigSq;
    }

    if( ( (M+1) % 10000 ) == 0){
      newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      Rcpp::Rcout << "iteration: " << M+1 << ": " << ctime(&newt) << "\n";
      Rcpp::checkUserInterrupt(); //checks if the user hit the "stop" icon to cancel running sampler.
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("samples") = Rcpp::List::create(
      Rcpp::Named("beta") = sample_beta.t(),
      Rcpp::Named("beta_tv") = sample_btv.t(),
      Rcpp::Named("mu") = sample_mu,
      Rcpp::Named("sigSq") = sample_sigSq),
    Rcpp::Named("accept") = Rcpp::List::create(
      Rcpp::Named("beta") = accept_beta,
      Rcpp::Named("beta_tv") = accept_btv,
      Rcpp::Named("mu") = accept_mu,
      Rcpp::Named("sigSq") = accept_sigSq)
  );

}

