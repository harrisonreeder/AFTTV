#ifndef AFT_LN_updates_H
#define AFT_LN_updates_H

double AFT_LN_loglik(const arma::mat &Wmat,
                     const arma::vec &wUInf,
                     const arma::vec &wLUeq,
                     const arma::vec &c0Inf,
                     const arma::mat &Xmat,
                     arma::vec &beta,
                     double &mu,
                     double &sigSq);

void AFT_LN_update_beta(const arma::mat &Wmat,
                 const arma::vec &wUInf,
                 const arma::vec &wLUeq,
                 const arma::vec &c0Inf,
                 const arma::mat &Xmat,
                 arma::vec &beta,
                 double &mu,
                 double &sigSq,
                 double &beta_prop_var,
                 arma::uvec &accept_beta,
                 double &curr_loglik);

void AFT_LN_update_mu(const arma::mat &Wmat,
               const arma::vec &wUInf,
               const arma::vec &wLUeq,
               const arma::vec &c0Inf,
               const arma::mat &Xmat,
               arma::vec &beta,
               double &mu,
               double &sigSq,
               double &mu_prop_var,
               int &accept_mu,
               double &curr_loglik);

void AFT_LN_update_sigSq(const arma::mat &Wmat,
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
               int &accept_sigSq,
               double &curr_loglik);

#endif

