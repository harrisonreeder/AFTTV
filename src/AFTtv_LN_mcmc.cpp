
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "AFTtv_utilities.h"
#include "AFTtv_LN_updates.h"
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
Rcpp::List AFTtv_LN_mcmc(const arma::mat &Ymat,
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
  int move; //index for which parameter to update
  int M; //counter for MCMC sampler
  int StoreInx; //index for where to store a sample, post-thinning


  //SET PRIOR CODES
  int sigSq_prior = prior_vec_num(1);
  int btv_prior = prior_vec_num(3);
  // int sigSq_prior = prior_list_num["sigSq"];
  // int btv_prior = prior_list_num["beta_tv"];

  //Rcpp::Rcout << "finished setting prior codes" << "\n";


  //SET MCMC TUNING PARAMETERS
  double mu_prop_var = tuning_vec(0); //might throw an error
  double sigSq_prop_var = tuning_vec(1); //might throw an error
  double beta_prop_var = tuning_vec(2); //might throw an error
  double btv_prop_var = tuning_vec(3); //might throw an error
  int K_max = tuning_vec(4); //maximum number of distinct piecewise intervals

  // double mu_prop_var = tuning_list["mu"]; //might throw an error
  // double sigSq_prop_var = tuning_list["sigSq"]; //might throw an error
  // double beta_prop_var = tuning_list["beta"]; //might throw an error
  // Rcpp::NumericVector btv_tuning = tuning_list["beta_tv"];
  // double btv_prop_var = btv_tuning[0]; //might throw an error
  // int K_max = btv_tuning[1]; //maximum number of distinct piecewise intervals

  //Rcpp::Rcout << "finished setting MCMC tuning params" << "\n";

  //SET HYPERPARAMETERS
  //flat prior for lognormal mu
  //flat prior for betas

  //gamma prior for lognormal sigSq
  double a_sigSq, b_sigSq;
  if(sigSq_prior==1){
    a_sigSq = hyper_vec(0);
    a_sigSq = hyper_vec(1);
    // Rcpp::NumericVector sigSq_hyper = hyper_list["sigSq"];
    // a_sigSq = sigSq_hyper[0];
    // b_sigSq = sigSq_hyper[1];
  }

  //mvn-icar prior for btv
  double a_btv, b_btv, c_btv,meanbtv,varbtv;
  arma::mat IminusW(K_max,K_max,arma::fill::eye);
  arma::vec Qvec(K_max);
  arma::mat Sigma_btv(K_max,K_max);
  arma::mat invSigma_btv(K_max,K_max);
  arma::mat cholinvSigma_btv(K_max,K_max);
  if(btv_prior==2){
    //flat hyperprior on mean for mean_btv
    //gamma hyperprior for var_btv
    a_btv = hyper_vec(2); //(This assumes that there is also an inv-gamma prior on sigSq, indexwise)
    b_btv = hyper_vec(3);
    // Rcpp::NumericVector btv_hyper = hyper_list["beta_tv"];
    // a_btv = btv_hyper[0];
    // b_btv = btv_hyper[1];
    //hyperparameter for dependence of time-varying beta coefficients
    c_btv = 1;

    //initialize matrices used in MVN-ICAR specification
    update_icar_mats(invSigma_btv, cholinvSigma_btv, IminusW, Qvec, knots_init, c_btv, K);

    //initialize hyperparameters
    meanbtv = start_vec[2+p+K];
    varbtv = start_vec[3+p+K];
    // meanbtv = start_list["meanbtv"];
    // varbtv = start_list["varbtv"];
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
  // double mu = start_list["mu"];
  // double sigSq = start_list["sigSq"];
  // arma::vec beta = start_list["beta"]; //TODO: what if there are no betas
  // arma::vec beta_tv = start_list["beta_tv"];

  arma::vec knots = knots_init;

  /* //next, we'll add these changing dimensional things:
   arma::vec beta_tv(K_max);
   beta_tv(span(0,K)) = startValues(arma::span(2+p,2+p-1+K));
   arma::vec knots(K_max);
   arma::vec knots(span(0,K)) = knots_init;
   */

  //Rcpp::Rcout << "finished setting starting values" << "\n";


  //create storage objects for results
  arma::mat sample_beta = arma::mat(p,n_store); //store betas column by column for "speed" ?
  arma::mat sample_btv = arma::mat(K,n_store); //store betas column by column for "speed" ?
  arma::vec sample_mu = arma::vec(n_store);
  arma::vec sample_sigSq = arma::vec(n_store);
  arma::vec accept_beta = arma::vec(p);
  arma::vec accept_btv = arma::vec(K);
  //initialize storage container for hyperparameters even if they aren't used.
  arma::vec sample_meanbtv = arma::vec(n_store); //gibbs sampled so no acceptance counter
  arma::vec sample_varbtv = arma::vec(n_store); //gibbs sampled so no acceptance counter
  int accept_mu = 0;
  int accept_sigSq = 0;

  //initialize V(t) values for yL, yU and cL
  arma::vec logVyL = Vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  arma::vec logVyU = Vx_pw(Ymat.col(1), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  arma::vec logVc0 = Vx_pw(Ymat.col(2), Xmat, beta, Xvec_tv, beta_tv, knots, 1);
  //initialize v(t) values for yL (only used when observation is exactly observed)
  arma::vec logvyL = vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv, knots, 1);

  //  Rcpp::Rcout << "initial logVyL: " << logVyL(arma::span(0,10)) << "\n";

  double curr_loglik = AFTtv_LN_loglik(logVyL,logVyU,logVc0,logvyL,
                                       yUInf,yLUeq,c0Inf,
                                       mu,sigSq);

  for(M = 0; M < n_iter; M++){
    //Rcpp::Rcout << "iteration: " << M << "\n";

    AFTtv_LN_update_mu(logVyL, logVyU, logVc0, logvyL,
                       yUInf, yLUeq, c0Inf,
                       mu, sigSq, mu_prop_var, accept_mu, curr_loglik);
    //Rcpp::Rcout << "updated mu: " << mu << "\n";
    AFTtv_LN_update_sigSq(logVyL, logVyU, logVc0, logvyL,
                          yUInf, yLUeq, c0Inf,
                          mu, sigSq, sigSq_prop_var,
                          a_sigSq, b_sigSq, accept_sigSq, curr_loglik);
    //Rcpp::Rcout << "updated sigSq: " << sigSq << "\n";
    AFTtv_LN_update_beta(logVyL, logVyU, logVc0, logvyL,
                         yUInf, yLUeq, c0Inf,
                         Xmat, beta, mu, sigSq, beta_prop_var, accept_beta, curr_loglik);
    //Rcpp::Rcout << "updated beta: " << beta.t() << "\n";


    if(btv_prior==2){ //use the ICAR sampler
      AFTtv_LN_update_hyper_icar(beta_tv, meanbtv, varbtv,
                                 a_btv, b_btv, cholinvSigma_btv);
      //Rcpp::Rcout << "updated hyper: mu" << meanbtv << "var" << varbtv << "\n";
      AFTtv_LN_update_btv_icar(Ymat, logVyL, logVyU, logVc0, logvyL,
                               yUInf, yLUeq, c0Inf,
                               Xmat, beta, Xvec_tv, beta_tv,
                               mu, sigSq, knots, btv_prop_var,
                               meanbtv, varbtv, cholinvSigma_btv, accept_btv, curr_loglik);
    } else{
      AFTtv_LN_update_btv(Ymat, logVyL, logVyU, logVc0, logvyL,
                               yUInf, yLUeq, c0Inf,
                               Xmat, beta, Xvec_tv, beta_tv,
                               mu, sigSq, knots, btv_prop_var, accept_btv, curr_loglik);
    }
    //Rcpp::Rcout << "updated beta_tv" << beta_tv.t() << "\n";

    // Storing posterior samples
    if( ( (M+1) % thin ) == 0 && (M+1) > n_burnin)
    {
      StoreInx = (M+1 - n_burnin)/thin;

      sample_beta.col(StoreInx - 1) = beta;
      sample_btv.col(StoreInx - 1) = beta_tv;
      sample_mu(StoreInx - 1) = mu;
      sample_sigSq(StoreInx - 1) = sigSq;
      if(btv_prior == 2){
        sample_meanbtv(StoreInx - 1) = meanbtv;
        sample_varbtv(StoreInx - 1) = varbtv;
      }
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
      Rcpp::Named("beta_tv") = sample_btv.t(),
      Rcpp::Named("meanbtv") = sample_meanbtv,
      Rcpp::Named("varbtv") = sample_varbtv),
      Rcpp::Named("accept") = Rcpp::List::create(
        Rcpp::Named("beta") = accept_beta,
        Rcpp::Named("beta_tv") = accept_btv,
        Rcpp::Named("mu") = accept_mu,
        Rcpp::Named("sigSq") = accept_sigSq)
  );

}

