
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
// [[Rcpp::depends(RcppArmadillo)]]


//COMPUTING V FUNCTIONS

double logsumexpnew(double x,double y){
  double c = fmax(x,y);
  return c + log(exp(x-c) + exp(y-c));
}

arma::vec v_pwnew(const arma::vec &T,
                  const arma::mat &Xmat,
                  const arma::vec &beta,
                  const arma::vec &Xvec_tv,
                  const arma::vec &beta_tv,
                  const arma::vec &knots,
                  const int log_out){
  // this function evaluates the v0 function, v0(t) = exp(-Xbeta(t))

  int n = T.n_rows;
  int K = knots.n_rows;
  arma::vec xbeta = Xmat * beta;
  arma::vec knots_newi = knots;
  arma::vec out(n);
  int i, k; //indices for nested for-loops
  int ind; //index

  for(i = 0; i < n; i++){ //for each observation, find the value of v at that time
    knots_newi(0) = exp(Xvec_tv(i)*beta_tv(0))*knots(0);
    for(k = 1; k < K; k++){
      knots_newi(k) = knots_newi(k-1) + exp(Xvec_tv(i)*beta_tv(k))*(knots(k)-knots(k-1));
    }
    knots_newi = knots_newi * exp(xbeta(i));

    ind = 0;
    while(T(i) >= knots_newi(ind)){ //should this be equal to or just greater than??
      ind++;
    }
    out(i) = - Xvec_tv(i) * beta_tv(ind); //initialize to v of interval
  }
  out = out - xbeta; // subtract time-invariant linear predictor

  if(log_out == 0){
    out = arma::exp(out);
  }
  return out;
}

arma::vec V_pwnew(const arma::vec &T,
                   const arma::mat &Xmat,
                   const arma::vec &beta,
                   const arma::vec &Xvec_tv,
                   const arma::vec &beta_tv,
                   const arma::vec &knots,
                   const int log_out){
  // this function evaluates the V0 function, V0(t) = int_0^t exp(-Xbeta(u)) du
  //and then multiplies it by the the baseline linear predictor

  int n = T.n_rows;
  int K = knots.n_rows; //check this!! and come up with rules for representation of knots
  arma::vec xbeta = Xmat * beta;
  arma::vec knots_newi = knots;
  arma::vec out(n);
  int i,k; //indices for nested for-loops
  double Del;

  for(i = 0; i < n; i++){ //for each observation
    // Rcpp::Rcout << "iteration: " << i << ": " << "\n";

    //if it's 0 (like for left-truncation time),
    //then the transformation will also be 0 and it also doesn't matter
    if(T(i)==0){
      continue;
    }

    // Rcpp::Rcout << "new knots: " << knots_newi << "\n";
    knots_newi(0) = exp(Xvec_tv(i)*beta_tv(0))*knots(0);
    for(k = 1; k < K; k++){
      // Rcpp::Rcout << "k: " << k << "\n";
      // Rcpp::Rcout << "new knots: " << knots_newi << "\n";
      knots_newi(k) = knots_newi(k-1) + exp(Xvec_tv(i)*beta_tv(k))*(knots(k)-knots(k-1));
    }
    knots_newi = knots_newi * exp(xbeta(i));
    // Rcpp::Rcout << "new knots: " << knots_newi << "\n";


    //loop through all of the intervals...
    //I exclude the 0 boundary knot from the 'knots' vector,
    //so kick things off with the length of the interval 0 to k_1 because we know it should be nonzero
    Del = fmax(0, fmin(knots_newi(0),T(i)));
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
      Del = fmax(0, fmin(knots_newi(k),T(i)) - knots_newi(k-1));

      if(Del > 0){ //if the duration is nonzero, then time is still accruing

        if(log_out == 1){
          //This parameterization sets each beta to be the height of its own interval
          out(i) = logsumexpnew(out(i), log(Del) - Xvec_tv(i)*(beta_tv(k)));
          //This parameterization yields beta "relative to the first interval" as in the frequentist code I wrote
          //out(i) = logsumexpnew(out(i), log(Del) - Xvec_tv(i)*(beta_tv(k) + beta_tv(0)));
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
  if(log_out == 1){
    out = out - xbeta; // multiply by time-invariant linear predictor
  } else{
    out = out % arma::exp(-xbeta); // multiply by time-invariant linear predictor
  }
  return out;

}

