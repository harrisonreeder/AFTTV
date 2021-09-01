#ifndef AFTtv_utilities_H
#define AFTtv_utilities_H


void inplace_tri_mat_mult(const arma::mat &trimat,
                          arma::vec &x);

double dmvnrm_arma(const arma::vec &x,
                   double &mean,
                   double &var,
                   const arma::mat &cholinvSigma_sub,
                   int logd = 1);

void update_icar_mats(arma::mat &invSigma_btv,
                      arma::mat &cholinvSigma_btv,
                      arma::mat &IminusW,
                      arma::vec &Qvec,
                      const arma::vec &knots,
                      double &c_btv,
                      int K);

arma::vec vx_pw(const arma::vec &T,
                const arma::mat &Xmat,
                const arma::vec &beta,
                const arma::vec &Xvec_tv,
                const arma::vec &beta_tv,
                const arma::vec &knots,
                const int log_out);

arma::vec Vx_pw(const arma::vec &T,
                const arma::mat &Xmat,
                const arma::vec &beta,
                const arma::vec &Xvec_tv,
                const arma::vec &beta_tv,
                const arma::vec &knots,
                const int log_out);

arma::vec Vx_log1p(const arma::vec &T,
                   const arma::mat &Xmat,
                   const arma::vec &beta,
                   const arma::vec &Xvec_tv,
                   double &beta_tv);

#endif

