#ifndef AFTtvnew_LN_updates_H
#define AFTtvnew_LN_updates_H

double AFTtvnew_LN_loglik(const arma::vec &logVyL,
                       const arma::vec &logVyU,
                       const arma::vec &logVc0,
                       const arma::vec &logvyL,
                       const arma::vec &yUInf,
                       const arma::vec &yLUeq,
                       const arma::vec &c0Inf,
                       double &mu,
                       double &sigSq);

void AFTtvnew_LN_update_mu(const arma::vec &logVyL,
                        const arma::vec &logVyU,
                        const arma::vec &logVc0,
                        const arma::vec &logvyL,
                        const arma::vec &yUInf,
                        const arma::vec &yLUeq,
                        const arma::vec &c0Inf,
                        double &mu,
                        double &sigSq,
                        double &mu_prop_var,
                        int &accept_mu,
                        double &curr_loglik);

void AFTtvnew_LN_update_sigSq(const arma::vec &logVyL,
                           const arma::vec &logVyU,
                           const arma::vec &logVc0,
                           const arma::vec &logvyL,
                           const arma::vec &yUInf,
                           const arma::vec &yLUeq,
                           const arma::vec &c0Inf,
                           double &mu,
                           double &sigSq,
                           double &sigSq_prop_var,
                           double &a_sigSq,
                           double &b_sigSq,
                           int &accept_sigSq,
                           double &curr_loglik);

void AFTtvnew_LN_update_beta(const arma::mat &Ymat,
                          arma::vec &VyL,
                          arma::vec &VyU,
                          arma::vec &Vc0,
                          arma::vec &vyL,
                          const arma::vec &yUInf,
                          const arma::vec &yLUeq,
                          const arma::vec &c0Inf,
                          const arma::mat &Xmat,
                          arma::vec &beta,
                          const arma::mat &Xvec_tv,
                          const arma::vec &beta_tv,
                          double &mu,
                          double &sigSq,
                          const arma::vec &knots,
                          double &beta_prop_var,
                          arma::vec &accept_beta,
                          double &curr_loglik);

void AFTtvnew_LN_update_btv(const arma::mat &Ymat,
                         arma::vec &VyL,
                         arma::vec &VyU,
                         arma::vec &Vc0,
                         arma::vec &vyL,
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
                         arma::vec &accept_btv,
                         double &curr_loglik);

#endif

