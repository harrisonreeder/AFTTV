
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
// [[Rcpp::depends(RcppArmadillo)]]

double logsumexp(double x,double y){
  double c = fmax(x,y);
  return c + log(exp(x-c) + exp(y-c));
}


arma::vec V0x_pw(const arma::vec &T,
                const arma::vec &Xvec_tv,
                const arma::vec &beta_tv,
                const arma::vec &knots,
                const int log_out){

  int n = T.n_rows;
  int K = knots.n_rows; //check this!! and come up with rules for representation of knots
  int i,k; //indices for nested for-loops
  double Del;
  arma::vec out(n);

  for(i = 0; i < n; i++){ //for each observation
    //Rcpp::Rcout << "iteration: " << i << ": " << "\n";

    //if it's 0 (like for left-truncation time),
    //then the transformation will also be 0 and it also doesn't matter
    if(T(i)==0){
      continue;
    }

    //loop through all of the intervals...
    //I exclude the 0 boundary knot from the 'knots' vector,
    //so kick things off with the length of the interval 0 to k_1 because we know it should be nonzero
    Del = fmax(0, fmin(knots(0),T(i)));
    //multiply the interval length by x_i*beta_tv
    if(log_out==1){
      out(i) = log(Del) - Xvec_tv(i)*beta_tv(0);
    } else{
      out(i) = Del*exp(-Xvec_tv(i)*beta_tv(0));
    }

    //this loop will only kick in if there are at least two non-zero knots,
    //meaning one besides the 'infinity' that ends the list of knots
    for(k = 1; k < K; k++){ //for each subsequent knot (k2, ..., kK)
      //duration of time spent in "k+1"th interval ("Del" for Delta as in SCR)
      Del = fmax(0, fmin(knots(k),T(i)) - knots(k-1));

      if(Del > 0){ //if the duration is nonzero, then time is still accruing

        if(log_out == 1){
          //This parameterization sets each beta to be the height of its own interval
          out(i) = logsumexp(out(i), log(Del) - Xvec_tv(i)*(beta_tv(k)));
          //This parameterization yields beta "relative to the first interval" as in the frequentist code I wrote
          //out(i) = logsumexp(out(i), log(Del) - Xvec_tv(i)*(beta_tv(k) + beta_tv(0)));
        } else{
          //This parameterization sets each beta to be the height of its own interval
          out(i) += Del*exp(-Xvec_tv(i)*beta_tv(k)); //add on kth element of x_i*beta_tv
          //This parameterization yields beta "relative to the first interval" as in the frequentist code I wrote
          //out(i) += Del*exp(-Xvec_tv(i)*(beta_tv(k)+beta_tv(0))); //add on kth element of x_i*beta_tv
        }

      } else{ //if the duration is zero, then time is up for this subject so move on
        break;
      }
    }
  }

  return out;

}

arma::vec v0x_pw(const arma::vec &T,
                 const arma::vec &Xvec_tv,
                 const arma::vec &beta_tv,
                 const arma::vec &knots,
                 const int log_out){
  int n = T.n_rows;
  int K = knots.n_rows; //check this!! and come up with rules for representation of knots
  int i,k; //indices for nested for-loops
  int ind; //index
  arma::vec out(n);

  for(i = 0; i < n; i++){ //for each observation, find the value of v at that time
    ind = 0;
    while(T(i) >= knots(ind)){ //should this be equal to or just greater than??
      ind++;
    }
    out(i) = - Xvec_tv(i) * beta_tv(ind); //initialize to v of first interval
  }
  if(log_out == 0){
    out = arma::exp(out);
  }
  return out;
}

// [[Rcpp::export]]
arma::vec Vx_pw(const arma::vec &T,
        const arma::mat &Xmat,
        const arma::vec &beta,
        const arma::vec &Xvec_tv,
        const arma::vec &beta_tv,
        const arma::vec &knots,
        const int log_out){
  arma::vec out = V0x_pw(T,Xvec_tv,beta_tv,knots,log_out);
  arma::vec xbeta = Xmat * beta;
  if(log_out == 1){
    out = out - xbeta; // multiply by time-invariant linear predictor
  } else{
    out = out % arma::exp(-xbeta); // multiply by time-invariant linear predictor
  }
  return out;
}

void AFTtv_LN_update_btv(const arma::mat &Ymat,
                          arma::vec &logVyL,
                          arma::vec &logVyU,
                          arma::vec &logVc0,
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
                          arma::vec &accept_btv){

  int n = logVyL.n_rows;
  int p_tv = beta_tv.n_rows;
  int j, i, u;
  double loglh, loglh_prop, logR;
  double logprior, logprior_prop;

  loglh = 0;
  loglh_prop = 0;

  //choose one beta_tv element to update
  //for simplicity, we will skip the ICAR and just go with a flat prior on these for now.
  //proposal is just a normal centered at previous value for now, nothing fancy there either.
  j = (int) R::runif(0, p_tv);
  arma::vec beta_tv_prop = beta_tv;
  beta_tv_prop(j) = R::rnorm(beta_tv(j), sqrt(btv_prop_var));

  //compute vectors with proposal V(t) values
  //later I could make this more efficient by subtracting and adding just the part that changed
  //but for now just recompute the whole thing
  arma::vec logVyL_prop = Vx_pw(Ymat.col(0), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);
  arma::vec logVyU_prop = Vx_pw(Ymat.col(1), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);
  arma::vec logVc0_prop = Vx_pw(Ymat.col(2), Xmat, beta, Xvec_tv, beta_tv_prop, knots,1);

  /* my problem is that the parameters for later intervals are flying off into space
  if(j>0){
    Rcpp::Rcout << "j : " << j ;
    Rcpp::Rcout << "beta_tvj_prop: " << beta_tv_prop(j) << "\n";
    Rcpp::Rcout << "initial logVyL: " << logVyL(arma::span(0,5)).t() << "\n";
    Rcpp::Rcout << "initial logVyL_prop: " << logVyL_prop(arma::span(0,5)).t() << "\n";
  }
  */

  for(i = 0; i < n; i++){
    if(yLUeq(i)==1){ //no censoring, so use log density
      loglh      += arma::log_normpdf(logVyL(i), mu, sqrt(sigSq));
      loglh_prop += arma::log_normpdf(logVyL_prop(i), mu, sqrt(sigSq));
    } else if(yUInf(i)==1) { //right censoring, so use log survival
      loglh      += log1p(-arma::normcdf(logVyL(i), mu, sqrt(sigSq)));
      loglh_prop += log1p(-arma::normcdf(logVyL_prop(i), mu, sqrt(sigSq)));
    } else { //interval censoring, so use log of difference of survival
      loglh      += log( arma::normcdf(logVyU(i), mu, sqrt(sigSq))
                           - arma::normcdf(logVyL(i), mu, sqrt(sigSq)) );
      loglh_prop += log( arma::normcdf(logVyU_prop(i), mu, sqrt(sigSq))
                           - arma::normcdf(logVyL_prop(i), mu, sqrt(sigSq)) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){
      loglh += -log1p(-arma::normcdf(logVc0(i), mu, sqrt(sigSq)));
      loglh_prop += -log1p(-arma::normcdf(logVc0_prop(i), mu, sqrt(sigSq)));
    }
  }
  logR = loglh_prop - loglh;

  //if we add a prior of some kind, it should go here...
  // here is a basic random walk prior, where the prior is centered at the previous value
  // I've just picked a default scale of 1 for the moment
   if(j>0){
   logprior = arma::log_normpdf(beta_tv(j), beta_tv(j-1), 1.0);
   logprior_prop = arma::log_normpdf(beta_tv_prop(j), beta_tv(j-1), 1.0);
   logR += logprior_prop - logprior;
   }


  //Rcpp::Rcout << "beta_tvj acceptance ratio: " << exp(logR) << "\n";

  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    beta_tv = beta_tv_prop;
    logVyL = logVyL_prop;
    logVyU = logVyU_prop;
    logVc0 = logVc0_prop;
    accept_btv(j) += 1;
  }
  return;
}

void AFTtv_LN_update_beta(arma::vec &logVyL,
                          arma::vec &logVyU,
                          arma::vec &logVc0,
                          const arma::vec &yUInf,
                          const arma::vec &yLUeq,
                          const arma::vec &c0Inf,
                          const arma::mat &Xmat,
                          arma::vec &beta,
                          double &mu,
                          double &sigSq,
                          double &beta_prop_var,
                          arma::vec &accept_beta){
  int n = logVyL.n_rows;
  int p = Xmat.n_cols;
  int j, i, u;
  double loglh, loglh_prop, logR;
  //double logprior, logprior_prop;

  loglh = 0;
  loglh_prop = 0;

  //update a single beta parameter
  j = (int) R::runif(0, p); //will need to update to sample only from non-tv X's
  double betaj_prop = R::rnorm(beta(j), sqrt(beta_prop_var));
  //compute updated V(t) values from the proposed betaj directly
  arma::vec xbetaj_propdiff = Xmat.col(j) * (beta(j) - betaj_prop);
  arma::vec logVyL_prop = logVyL + xbetaj_propdiff;
  arma::vec logVyU_prop = logVyU + xbetaj_propdiff;
  arma::vec logVc0_prop = logVc0 + xbetaj_propdiff;

  for(i = 0; i < n; i++){
    if(yLUeq(i)==1){ //no censoring, so use log density
      loglh      += arma::log_normpdf(logVyL(i), mu, sqrt(sigSq));
      loglh_prop += arma::log_normpdf(logVyL_prop(i), mu, sqrt(sigSq));
    } else if(yUInf(i)==1) { //right censoring, so use log survival
      loglh      += log1p(-arma::normcdf(logVyL(i), mu, sqrt(sigSq)));
      loglh_prop += log1p(-arma::normcdf(logVyL_prop(i), mu, sqrt(sigSq)));
    } else { //interval censoring, so use log of difference of survival
      loglh      += log( arma::normcdf(logVyU(i), mu, sqrt(sigSq))
                           - arma::normcdf(logVyL(i), mu, sqrt(sigSq)) );
      loglh_prop += log( arma::normcdf(logVyU_prop(i), mu, sqrt(sigSq))
                           - arma::normcdf(logVyL_prop(i), mu, sqrt(sigSq)) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){
      loglh += -log1p(-arma::normcdf(logVc0(i), mu, sqrt(sigSq)));
      loglh_prop += -log1p(-arma::normcdf(logVc0_prop(i), mu, sqrt(sigSq)));
    }
  }
  logR = loglh_prop - loglh;

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
    accept_beta(j) += 1;
  }
  return;
}

void AFTtv_LN_update_mu(const arma::vec &logVyL,
                        const arma::vec &logVyU,
                        const arma::vec &logVc0,
                        const arma::vec &yUInf,
                        const arma::vec &yLUeq,
                        const arma::vec &c0Inf,
                        double &mu,
                        double &sigSq,
                        double &mu_prop_var,
                        int &accept_mu) {
  int n = logVyL.n_rows;
  int i, u;
  double loglh, loglh_prop, logR, mu_prop;
  //double logprior, logprior_prop;

  loglh = 0;
  loglh_prop = 0;
  mu_prop = R::rnorm(mu, sqrt(mu_prop_var));

  //in the future, I'd love to turn the yUInf etc indicators into indices so I could
  //subset vectors and use vector operations directly.

  for(i = 0; i < n; i++){
    /* cases:
     * interval censoring (left and right diff.): compute survivor function at both, ll is log of difference.
     * (note: log of difference is also logsumexp of log survivor functions, in case it's useful)
     * right censoring (right is infinite): compute survivor function at left, ll is that
     * no censoring (left and right equal): compute log density, ll is that */
    if(yLUeq(i)==1){ //no censoring, so use log density
      loglh      += arma::log_normpdf(logVyL(i), mu, sqrt(sigSq));
      loglh_prop += arma::log_normpdf(logVyL(i), mu_prop, sqrt(sigSq));
    } else if(yUInf(i)==1) { //right censoring, so use log survival
      loglh      += log1p(-arma::normcdf(logVyL(i), mu, sqrt(sigSq)));
      loglh_prop += log1p(-arma::normcdf(logVyL(i), mu_prop, sqrt(sigSq)));
    } else { //interval censoring, so use log of difference of survival
      loglh      += log( arma::normcdf(logVyU(i), mu, sqrt(sigSq))
                           - arma::normcdf(logVyL(i), mu, sqrt(sigSq)) );
      loglh_prop += log( arma::normcdf(logVyU(i), mu_prop, sqrt(sigSq))
                           - arma::normcdf(logVyL(i), mu_prop, sqrt(sigSq)) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){
      loglh += -log1p(-arma::normcdf(logVc0(i), mu, sqrt(sigSq)));
      loglh_prop += -log1p(-arma::normcdf(logVc0(i), mu_prop, sqrt(sigSq)));
    }
  }

  logR = loglh_prop - loglh;

  /* //if we add a prior of some kind, it should go here...
   logprior = R::dnorm(mu, mu0, sqrt(h0), 1);
   logprior_prop = R::dnorm(mu_prop, mu0, sqrt(h0), 1);
   logR += logprior_prop - logprior;
   */

  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    mu = mu_prop;
    accept_mu += 1;
  }
  return;

}

void AFTtv_LN_update_sigSq(const arma::vec &logVyL,
                           const arma::vec &logVyU,
                           const arma::vec &logVc0,
                           const arma::vec &yUInf,
                           const arma::vec &yLUeq,
                           const arma::vec &c0Inf,
                           double &mu,
                           double &sigSq,
                           double &sigSq_prop_var,
                           double &a_sigSq,
                           double &b_sigSq,
                           int &accept_sigSq) {
  int n = logVyL.n_rows;
  int i, u;
  double loglh, loglh_prop, logR, sigSq_prop;
  double logprior, logprior_prop;

  loglh = 0;
  loglh_prop = 0;
  sigSq_prop = exp( R::rnorm( log(sigSq), sqrt(sigSq_prop_var) ) );

  for(i = 0; i < n; i++){
    if(yLUeq(i)==1){ //no censoring, so use log density
      loglh      += arma::log_normpdf(logVyL(i), mu, sqrt(sigSq));
      loglh_prop += arma::log_normpdf(logVyL(i), mu, sqrt(sigSq_prop));
    } else if(yUInf(i)==1) { //right censoring, so use log survival
      loglh      += log1p(-arma::normcdf(logVyL(i), mu, sqrt(sigSq)));
      loglh_prop += log1p(-arma::normcdf(logVyL(i), mu, sqrt(sigSq_prop)));
    } else { //interval censoring, so use log of difference of survival
      loglh      += log( arma::normcdf(logVyU(i), mu, sqrt(sigSq))
                           - arma::normcdf(logVyL(i), mu, sqrt(sigSq)) );
      loglh_prop += log( arma::normcdf(logVyU(i), mu, sqrt(sigSq_prop))
                           - arma::normcdf(logVyL(i), mu, sqrt(sigSq_prop)) );
    }
    //if there is left-truncation, then subtract off the extra term
    if(c0Inf(i) == 0){
      loglh += -log1p(-arma::normcdf(logVc0(i), mu, sqrt(sigSq)));
      loglh_prop += -log1p(-arma::normcdf(logVc0(i), mu, sqrt(sigSq_prop)));
    }
  }

  //sigSq has an inverse-gamma prior, so we define them as follows
  logprior = (-a_sigSq-1)*log(sigSq) - b_sigSq/sigSq;
  logprior_prop = (-a_sigSq-1)*log(sigSq_prop) - b_sigSq/sigSq_prop;

  //Last "log(sigSq_prop) - log(sigSq)" is the "jacobian" because we sample on transformed scale
  //I think because we put prior on sigSq, but then sample log(sigSq), so we need to add jacobian

  logR = loglh_prop - loglh + logprior_prop - logprior + log(sigSq_prop) - log(sigSq);
  u = log(R::runif(0, 1)) < logR;
  if(u == 1){
    sigSq = sigSq_prop;
    accept_sigSq += 1;
  }
  return;

}






//SOMEDAY THIS WILL BE USEFUL


arma::vec Vx_log1p(const arma::vec &T,
                   const arma::mat &Xmat,
                   const arma::vec &beta,
                   const arma::vec &Xvec_tv,
                   double &beta_tv){
  int n = T.n_rows;
  int i; //for-loop index
  arma::vec Vy(n);
  arma::vec xbeta = Xmat * beta;
  arma::vec xbeta_tv = Xvec_tv * beta_tv;

  for(i = 0; i < n; i++){ //for each observation
    if(xbeta_tv(i) == 1){
      Vy(i) = log1p(T(i));
    } else{
      Vy(i) = (pow(T(i)+1,1-xbeta_tv(i))-1)/(1-xbeta_tv(i));
    }
    Vy(i) = Vy(i) * exp(-xbeta(i)); // multiply by time-invariant linear predictor
  }
  return Vy;
}
