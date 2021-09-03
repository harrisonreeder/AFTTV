
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "AFT_LN_updates.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List AFT_LN_mcmc(const arma::mat& Wmat,
                        const arma::vec& wUInf,
                        const arma::vec& wLUeq,
                        const arma::vec& c0Inf,
                        const arma::mat& Xmat,
                        const arma::vec& hyper_vec,
                        const arma::vec& tuning_vec,
                        const arma::vec& start_vec,
                        int n_burnin,
                        int n_sample,
                        int thin){
  //timekeeping objects
  std::time_t newt;

  //set constants
  int p = Xmat.n_cols;
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'
  int move; //index for which parameter to update
  int M; //counter for MCMC sampler
  int StoreInx; //index for where to store a sample, post-thinning

  //set hyper_vecarameters
  double a_sigSq = hyper_vec(0);
  double b_sigSq = hyper_vec(1);

  //set MCMC tuning parameters
  double mu_prop_var = tuning_vec(0);
  double sigSq_prop_var = tuning_vec(1);
  double beta_prop_var = tuning_vec(2);

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
  double mu = start_vec(0);
  double sigSq = start_vec(1);
  arma::vec beta = start_vec(arma::span(2,2+p-1));

  //create storage objects
  arma::mat sample_beta = arma::mat(p,n_store); //store betas column by column for "speed" ?
  arma::vec sample_mu = arma::vec(n_store);
  arma::vec sample_sigSq = arma::vec(n_store);
  arma::uvec accept_beta = arma::uvec(p, arma::fill::zeros);
  int accept_mu = 0;
  int accept_sigSq = 0;

  double curr_loglik = AFT_LN_loglik(Wmat,wUInf,wLUeq,c0Inf,
                             Xmat,beta,mu,sigSq);

  for(M = 0; M < n_iter; M++){

    //if we've changed beta, recompute V(t) (or in time-invariant case just eta), otherwise no need to

    AFT_LN_update_beta(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, beta_prop_var, accept_beta,curr_loglik);
    AFT_LN_update_mu(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, mu_prop_var, accept_mu,curr_loglik);
    AFT_LN_update_sigSq(Wmat, wUInf, wLUeq, c0Inf, Xmat, beta, mu, sigSq, sigSq_prop_var, a_sigSq, b_sigSq, accept_sigSq,curr_loglik);

    /*
    move = (int) R::runif(0, 3); //for now, equal probability of each of the 3 moves.

    if(move == 0){
      AFT_LN_update_beta(Wmat, wUInf, wLUeq, c0Inf,
                         Xmat, beta, mu, sigSq,
                         beta_prop_var, accept_beta, curr_loglik);
    }
    if(move == 1){
      AFT_LN_update_mu(Wmat, wUInf, wLUeq, c0Inf,
                       Xmat, beta, mu, sigSq,
                       mu_prop_var, accept_mu, curr_loglik);
    }
    if(move == 2){
      AFT_LN_update_sigSq(Wmat, wUInf, wLUeq, c0Inf,
                          Xmat, beta, mu, sigSq,
                          sigSq_prop_var, a_sigSq, b_sigSq, accept_sigSq, curr_loglik);
    }
    */

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
      Rcpp::Named("mu") = sample_mu,
      Rcpp::Named("sigSq") = sample_sigSq,
      Rcpp::Named("beta") = sample_beta.t()),
    Rcpp::Named("accept") = Rcpp::List::create(
      Rcpp::Named("mu") = accept_mu,
      Rcpp::Named("sigSq") = accept_sigSq,
      Rcpp::Named("beta") = accept_beta)
  );

}

