
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "AFTtvnew_utilities.h"
#include "AFTtvnew_LN_updates.h"
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
Rcpp::List AFTtvnew_LN_mcmc(const arma::mat &Ymat,
                         const arma::vec &yUInf,
                         const arma::vec &yLUeq,
                         const arma::vec &c0Inf,
                         const arma::mat &Xmat,
                         const arma::vec &Xvec_tv,
                         const arma::vec &prior_vec_num,
                         const arma::vec &hyper_vec,
                         const arma::vec &tuning_vec,
                         const arma::vec &start_vec,
                         const arma::vec &knots_init,
                         int n_burnin,
                         int n_sample,
                         int thin){

  // Rcpp::List prior_list_num,
  // Rcpp::List hyper_list,
  // Rcpp::List tuning_list,
  // Rcpp::List start_list,


  //timekeeping objects
  std::time_t newt;

  //SET CONSTANTS
  int p = Xmat.n_cols;
  //one parameter for every element of the "knots" vector, because it excludes 0 but includes infty.
  int K = knots_init.n_rows;
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'
  // int move; //index for which parameter to update
  int M; //counter for MCMC sampler
  int StoreInx; //index for where to store a sample, post-thinning


  //SET PRIOR CODES
  int sigSq_prior = prior_vec_num(1);
  // int btv_prior = prior_vec_num(3);
  // int sigSq_prior = prior_list_num["sigSq"];
  // int btv_prior = prior_list_num["beta_tv"];

  //Rcpp::Rcout << "finished setting prior codes" << "\n";


  //SET MCMC TUNING PARAMETERS
  double mu_prop_var = tuning_vec(0); //might throw an error
  double sigSq_prop_var = tuning_vec(1); //might throw an error
  double beta_prop_var = tuning_vec(2); //might throw an error
  double btv_prop_var = tuning_vec(2); //might throw an error

  //Rcpp::Rcout << "finished setting MCMC tuning params" << "\n";

  //SET HYPERPARAMETERS
  //flat prior for lognormal mu
  //flat prior for betas

  //gamma prior for lognormal sigSq
  double a_sigSq, b_sigSq;
  if(sigSq_prior==1){
    a_sigSq = hyper_vec(0);
    a_sigSq = hyper_vec(1);
  }

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

  //Rcpp::Rcout << "finished setting hyperparams" << "\n";


  //initialize starting values
  double mu = start_vec(0);
  double sigSq = start_vec(1);
  arma::vec beta = start_vec(arma::span(2,1+p)); //TODO: what if there are no betas
  arma::vec beta_tv = start_vec(arma::span(2+p,1+p+K));

  arma::vec knots = knots_init;

  //Rcpp::Rcout << "finished setting starting values" << "\n";

  //create storage objects for results
  arma::mat sample_beta = arma::mat(p,n_store); //store betas column by column for "speed" ?
  arma::mat sample_btv = arma::mat(K,n_store); //store betas column by column for "speed" ?
  arma::vec sample_mu = arma::vec(n_store);
  arma::vec sample_sigSq = arma::vec(n_store);
  arma::vec accept_beta = arma::vec(p);
  arma::vec accept_btv = arma::vec(K);
  int accept_mu = 0;
  int accept_sigSq = 0;

  //initialize V(t) values for yL, yU and cL
  arma::vec logVyL = V_pwnew(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  arma::vec logVyU = V_pwnew(Ymat.col(1), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  arma::vec logVc0 = V_pwnew(Ymat.col(2), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  //initialize v(t) values for yL (only used when observation is exactly observed)
  arma::vec logvyL = v_pwnew(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv, knots, 1);

  //Rcpp::Rcout << "initial logVyL: " << logVyL(arma::span(0,10)) << "\n";
  //Rcpp::Rcout << "initial logvyL: " << logvyL(arma::span(0,10)) << "\n";

  double curr_loglik = AFTtvnew_LN_loglik(logVyL,logVyU,logVc0,logvyL,
                                       yUInf,yLUeq,c0Inf,
                                       mu,sigSq);

  for(M = 0; M < n_iter; M++){
    //Rcpp::Rcout << "iteration: " << M << "\n";

    AFTtvnew_LN_update_mu(logVyL, logVyU, logVc0, logvyL,
                       yUInf, yLUeq, c0Inf,
                       mu, sigSq, mu_prop_var, accept_mu, curr_loglik);
    //Rcpp::Rcout << "updated mu: " << mu << "\n";
    AFTtvnew_LN_update_sigSq(logVyL, logVyU, logVc0, logvyL,
                          yUInf, yLUeq, c0Inf,
                          mu, sigSq, sigSq_prop_var,
                          a_sigSq, b_sigSq, accept_sigSq, curr_loglik);
    //Rcpp::Rcout << "updated sigSq: " << sigSq << "\n";
    AFTtvnew_LN_update_beta(Ymat, logVyL, logVyU, logVc0, logvyL,
                         yUInf, yLUeq, c0Inf,
                         Xmat, beta, Xvec_tv, beta_tv,
                         mu, sigSq, knots, beta_prop_var, accept_beta, curr_loglik);
    //Rcpp::Rcout << "updated beta: " << beta.t() << "\n";

    AFTtvnew_LN_update_btv(Ymat, logVyL, logVyU, logVc0, logvyL,
                               yUInf, yLUeq, c0Inf,
                               Xmat, beta, Xvec_tv, beta_tv,
                               mu, sigSq, knots, btv_prop_var, accept_btv, curr_loglik);
    //Rcpp::Rcout << "updated beta_tv" << beta_tv.t() << "\n";

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
      Rcpp::Named("mu") = sample_mu,
      Rcpp::Named("sigSq") = sample_sigSq,
      Rcpp::Named("beta") = sample_beta.t(),
      Rcpp::Named("beta_tv") = sample_btv.t()),
      Rcpp::Named("accept") = Rcpp::List::create(
        Rcpp::Named("beta") = accept_beta,
        Rcpp::Named("beta_tv") = accept_btv,
        Rcpp::Named("mu") = accept_mu,
        Rcpp::Named("sigSq") = accept_sigSq)
  );

}

