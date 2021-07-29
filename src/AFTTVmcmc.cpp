
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List BAFTtvLTmcmc(const arma::mat& Wmat,
                        const arma::vec& wUInf,
                        const arma::vec& c0Inf,
                        const arma::mat& Xmat,
                        const arma::vec& hyperP,
                        const arma::vec& mcmcP,
                        const arma::vec& startValues,
                        int n_burnin,
                        int n_sample,
                        int thin){
  //set constants
  int n = Xmat.n_rows;
  int p = Xmat.n_cols;
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'
  int move; //index for which parameter to update
  int M; //counter for MCMC sampler

  //set hyperparameters
  double a_sigSq = hyperP(0);
  double b_sigSq = hyperP(1);

  //set MCMC tuning parameters
  double beta_prop_var = mcmcP(0);
  double mu_prop_var = mcmcP(1);
  double sigSq_prop_var = mcmcP(2);

  /*
  double pBeta = 0.3; //probability of updating beta (random scan)
  double pMu = 0.3; //probability of updating mu (random scan)
  double pSigSq = 0.2;  //probability of updating sigSq (random scan)

  arma::vec moveProb = arma::vec(6);
  moveProb(0) = pBeta;
  moveProb(1) = pMu;
  moveProb(2) = pSigSq;
  */


  //initialize starting values
  double mu = startValues(0);
  double sigSq = startValues(1);
  arma::vec beta = startValues(arma::span(2,2+p-1));

  //create storage objects
  arma::mat sample_beta = arma::mat(n_store,p);
  arma::vec sample_mu = arma::vec(n_store);
  arma::vec sample_sigSq = arma::vec(n_store);
  arma::vec accept_beta = arma::vec(p);
  double accept_mu = 0;
  double accept_sigSq = 0;

  for(M = 0; M < n_iter; M++){
    /*
    if( ( (M) % 100 ) == 0)
    {
      Rcpp::checkUserInterrupt();
    }
    */

    // move = R::sample(rr, moveProb, 5);

    move = 1111;

    //update_beta(c0, c0_neginf, X, w, beta, tauSq, mu, sigSq, beta_prop_var, accept_beta);
    //update_mu(c0, c0_neginf, X, w, beta, &mu, sigSq, mu0, h0, mu_prop_var, &accept_mu);
    //update_sigSq(c0, c0_neginf, X, w, beta, tauSq, mu, &sigSq, a_sigSq, b_sigSq, sigSq_prop_var, &accept_sigSq);

  }

  return Rcpp::List::create(
    Rcpp::Named("samples") = Rcpp::List::create(
      Rcpp::Named("beta") = sample_beta,
      Rcpp::Named("mu") = sample_mu,
      Rcpp::Named("sigSq") = sample_sigSq),
    Rcpp::Named("accept") = Rcpp::List::create(
      Rcpp::Named("beta") = accept_beta,
      Rcpp::Named("mu") = accept_mu,
      Rcpp::Named("sigSq") = accept_sigSq)
  );

}

