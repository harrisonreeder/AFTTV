
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
// [[Rcpp::depends(RcppArmadillo)]]



//MVN DENSITY (adapted from https://gallery.rcpp.org/articles/dmvnorm_arma/)
static double const log2pi = std::log(2.0 * arma::datum::pi);

void inplace_tri_mat_mult(const arma::mat &trimat,
                          arma::vec &x){
  //returns inplace R * x for x, where R is upper triangular
  int n = trimat.n_cols;

  for(int i = 0; i < n; i++){
    double tmp(0.);
    for(int j = i; j < n; j++){
      tmp += trimat.at(i,j) * x(j);
    }
    x(i) = tmp;
  }
  return;
}

// [[Rcpp::export]]
double dmvnrm_arma(const arma::vec &x,
                   double &mean,
                   double &var,
                   const arma::mat &cholinvSigma_sub,
                   int logd = 1) {
  //this assumes that the covariance matrix is var * Sigma, but you input var and invSigma
  //for now, x is a matrix where each of the n rows is MVN, each of dimension xdim.
  double xdim = 1.0 * x.n_cols;
  double out;
  //Now, get log of 1/determinant(Sigma) by the method of
  //https://math.stackexchange.com/questions/3158303/using-cholesky-decomposition-to-compute-covariance-matrix-determinant
  //finally, note that det(1/var * invSigma) = (1/var)^(K) det(invSigma) https://math.stackexchange.com/questions/1215555/determinant-of-matrix-times-a-constant
  //so, det(var*Sigma)^(-1/2) = var^(-K/2) det(Sigma)^(-1/2) = var^(-K/2) det(invSigma)^(1/2)
  // = -(K/2)*log(var) + logdet(invSigma)/2 = -(K/2)*log(var) + sum(log(chol(invSigma)))
  double other_terms = arma::sum(arma::log(cholinvSigma_sub.diag()))
                       - xdim/2.0 * (log2pi + log(var));

  arma::vec z;
  z = (x - mean);
  inplace_tri_mat_mult(cholinvSigma_sub, z);
  //now, <z,z> = (x-mean).t() * R.t() * R * (x-mean) = (x-mean).t() * invSigma * (x-mean)
  out = other_terms - 0.5 * arma::dot(z, z) / var;

  if (logd==1){
    return out;
  } else{
    return exp(out);
  }
}

//COMPUTING MVN-ICAR MATRICES

void update_icar_mats(arma::mat &invSigma_btv,
                      arma::mat &cholinvSigma_btv,
                      arma::mat &IminusW,
                      arma::vec &Qvec,
                      const arma::vec &knots,
                      double &c_btv,
                      int K){

  //Three cases that affect how to compute Q and W:
    //there's exactly 1 interval
    //there's exactly 2 intervals
    //there's 3 or more intervals

  //because there's no 0, then the first 'interval' length is just the first knot.
  //THIS WILL ONLY MAKE SENSE IF THE "ENDPOINT" IS ACTUALLY THE LARGEST OBSERVED VALUE AND NOT INFINITY
  arma::vec intervals(K);
  intervals(0) = knots(0);

  Rcpp::Rcout << "knots" << knots << "\n";


  if(K==1){
    Qvec(0) = 1/intervals(0); //"length" of interval is maximum observed time;
    IminusW(0,0) = 1; //no weights because no neighboring intervals
    //Sigma_btv(0,0) = Qvec(0);
    invSigma_btv(0,0) = 1/Qvec(0);
  } else{
    //fill in remaining intervals after the first, using diff function to get lengths
    //end at K-1 because there is no Kth element due to 0 indexing
    intervals(arma::span(1,K-1)) = arma::diff(knots(arma::span(0,K-1)));

    //Rcpp::Rcout << "intervals" << intervals << "\n";

    //fill in two endpoint Q values
    Qvec(0) =   2/(2*intervals(0) + intervals(1));
    Qvec(K-1) = 2/(intervals(K-2) + 2*intervals(K-1));
    IminusW(0,1) = -c_btv * (intervals(0) + intervals(1)) * Qvec(0) / 2;
    IminusW(K-1,K-2) = -c_btv * (intervals(K-1) + intervals(K-2)) * Qvec(K-1) / 2;
    //this loop will kick in if there are more than two intervals
    for(int k=1; k<(K-1); k++){
      Qvec(k) = 2/(intervals(k-1) + 2*intervals(k) + intervals(k+1));
      IminusW(k,k-1) = -c_btv * (intervals(k-1) + intervals(k)) * Qvec(k) / 2;
      IminusW(k,k+1) = -c_btv * (intervals(k) + intervals(k+1)) * Qvec(k) / 2;
    }

    //Rcpp::Rcout << "Here is IminusW " << IminusW(arma::span(0,K-1),arma::span(0,K-1)) << "\n";


    //Sigma = (I - W)^{-1} diag(Q)
    //arma::mat Sigma_btv = solve(IminusW(arma::span(0,K-1),arma::span(0,K-1)), arma:: diagmat(Qvec(arma::span(0,K-1))));
    //Rcpp::Rcout << "Here is Sigma_btv " << Sigma_btv << "\n";

    //Sigma_btv(arma::span(0,K-1),arma::span(0,K-1)) = solve(IminusW(arma::span(0,K-1),arma::span(0,K-1)),  diagmat(Qvec(arma::span(0,K-1))));
    //Sigma_btv(arma::span(0,K-1),arma::span(0,K-1)) = inv(IminusW(arma::span(0,K-1),arma::span(0,K-1))).each_row() % Qvec(arma::span(0,K-1));

    //Identity for square matrices (AB)^{-1} = B^{-1}A^{-1} (https://proofwiki.org/wiki/Inverse_of_Matrix_Product)
    //So, invSigma = diag(1/Q) (I-W)
    //This in turn means just multiplying every column of (I-W) by 1/Q
    invSigma_btv(arma::span(0,K-1),arma::span(0,K-1)) = IminusW(arma::span(0,K-1),arma::span(0,K-1)).each_col() % (1/Qvec(arma::span(0,K-1)));
    //invSigma_btv(arma::span(0,K-1),arma::span(0,K-1)) = arma::inv(Sigma_btv(arma::span(0,K-1),arma::span(0,K-1)));

    //Rcpp::Rcout << "Here is the invSigma_btv " << invSigma_btv(arma::span(0,K-1),arma::span(0,K-1)) << "\n";


    //Identity for square matrices (AB)^{-1} = B^{-1}A^{-1} (https://proofwiki.org/wiki/Inverse_of_Matrix_Product)
    //So, invSigma = arma::diag(1/Q) (I-W)
    //This in turn means just multiplying every column of (I-W) by 1/Q
    //Add on a cholesky decomposition and we get upper triangular R s.t. RtR = invSigma
    cholinvSigma_btv(arma::span(0,K-1),arma::span(0,K-1)) = arma::chol( invSigma_btv(arma::span(0,K-1),arma::span(0,K-1)) );
  }

  Rcpp::Rcout << "Here is the cholinvSigma_btv " << cholinvSigma_btv(arma::span(0,K-1),arma::span(0,K-1)) << "\n";

  return;

}



//COMPUTING V FUNCTIONS

double logsumexp(double x,double y){
  double c = fmax(x,y);
  return c + log(exp(x-c) + exp(y-c));
}

arma::vec v0x_pw(const arma::vec &T,
                 const arma::vec &Xvec_tv,
                 const arma::vec &beta_tv,
                 const arma::vec &knots,
                 const int log_out){
  // this function evaluates the v0 function, v0(t) = exp(-Xbeta(t))

  int n = T.n_rows;
  int i; //indices for nested for-loops
  int ind; //index
  arma::vec out(n);

  for(i = 0; i < n; i++){ //for each observation, find the value of v at that time
    ind = 0;
    while(T(i) >= knots(ind)){ //should this be equal to or just greater than??
      ind++;
    }
    out(i) = - Xvec_tv(i) * beta_tv(ind); //initialize to v of interval
  }
  if(log_out == 0){
    out = arma::exp(out);
  }
  return out;
}

arma::vec vx_pw(const arma::vec &T,
                const arma::mat &Xmat,
                const arma::vec &beta,
                const arma::vec &Xvec_tv,
                const arma::vec &beta_tv,
                const arma::vec &knots,
                const int log_out){
  // this function evaluates the v function, v0(t) = exp(Xbeta - Xbeta(t))

  arma::vec out = v0x_pw(T,Xvec_tv,beta_tv,knots,log_out);
  arma::vec xbeta = Xmat * beta;
  if(log_out == 1){
    out = out - xbeta; // subtract time-invariant linear predictor
  } else{
    out = out % arma::exp(-xbeta); // multiply exponentiated time-invariant linear predictor
  }
  return out;
}




arma::vec V0x_pw(const arma::vec &T,
                const arma::vec &Xvec_tv,
                const arma::vec &beta_tv,
                const arma::vec &knots,
                const int log_out){
  // this function evaluates the V0 function, V0(t) = int_0^t exp(-Xbeta(u)) du

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

// [[Rcpp::export]]
arma::vec Vx_pw(const arma::vec &T,
        const arma::mat &Xmat,
        const arma::vec &beta,
        const arma::vec &Xvec_tv,
        const arma::vec &beta_tv,
        const arma::vec &knots,
        const int log_out){
  // this function evaluates the V function, V(t) = int_0^t exp(-Xbeta - Xbeta(u)) du

  arma::vec out = V0x_pw(T,Xvec_tv,beta_tv,knots,log_out);
  arma::vec xbeta = Xmat * beta;
  if(log_out == 1){
    out = out - xbeta; // multiply by time-invariant linear predictor
  } else{
    out = out % arma::exp(-xbeta); // multiply by time-invariant linear predictor
  }
  return out;
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

