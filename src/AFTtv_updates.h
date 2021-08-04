#ifndef AFTtv_updates_H
#define AFTtv_updates_H

void update_beta(const arma::mat &Wmat,
                 const arma::vec &wUInf,
                 const arma::vec &wLUeq,
                 const arma::vec &c0Inf,
                 const arma::mat &Xmat,
                 arma::vec &beta,
                 double &mu,
                 double &sigSq,
                 double &beta_prop_var,
                 arma::vec &accept_beta);

void update_mu(const arma::mat &Wmat,
               const arma::vec &wUInf,
               const arma::vec &wLUeq,
               const arma::vec &c0Inf,
               const arma::mat &Xmat,
               arma::vec &beta,
               double &mu,
               double &sigSq,
               double &mu_prop_var,
               int &accept_mu);

void update_sigSq(const arma::mat &Wmat,
               const arma::vec &wUInf,
               const arma::vec &wLUeq,
               const arma::vec &c0Inf,
               const arma::mat &Xmat,
               arma::vec &beta,
               double &mu,
               double &sigSq,
               double &sigSq_prop_var,
               double &a_sigSq,
               double &b_sigSq,
               int &accept_sigSq);

#endif

