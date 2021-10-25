#ifndef AFTtvnew_utilities_H
#define AFTtvnew_utilities_H

arma::vec v_pwnew(const arma::vec &T,
                const arma::mat &Xmat,
                const arma::vec &beta,
                const arma::vec &Xvec_tv,
                const arma::vec &beta_tv,
                const arma::vec &knots,
                const int log_out);

arma::vec V_pwnew(const arma::vec &T,
                const arma::mat &Xmat,
                const arma::vec &beta,
                const arma::vec &Xvec_tv,
                const arma::vec &beta_tv,
                const arma::vec &knots,
                const int log_out);

#endif

