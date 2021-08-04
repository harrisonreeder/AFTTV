
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "AFTtv_LN_updates.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List AFTtv_LN_mcmc(const arma::mat& Wmat,
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
    AFTtv_LN_update_beta(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, beta_prop_var, accept_beta);
    AFTtv_LN_update_mu(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, mu_prop_var, accept_mu);
    AFTtv_LN_update_sigSq(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, sigSq_prop_var, a_sigSq, b_sigSq, accept_sigSq);
    */

    move = (int) R::runif(0, 3); //for now, equal probability of each of the 3 moves.

    if(move == 0){
      AFTtv_LN_update_beta(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, beta_prop_var, accept_beta);
    }
    if(move == 1){
      AFTtv_LN_update_mu(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, mu_prop_var, accept_mu);
    }
    if(move == 2){
      AFTtv_LN_update_sigSq(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, sigSq_prop_var, a_sigSq, b_sigSq, accept_sigSq);
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

