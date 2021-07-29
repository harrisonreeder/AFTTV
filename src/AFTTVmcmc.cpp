
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
  int n_store = n_sample / thin; //tested elsewhere that this 'fits'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'

  //set hyperparameters
  double a_sigSq = hyperP(0);
  double b_sigSq = hyperP(1);

  //set MCMC tuning parameters
  double beta_prop_var = mcmcP(0);
  double mu_prop_var = mcmcP(1);
  double sigSq_prop_var = mcmcP(2);

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

