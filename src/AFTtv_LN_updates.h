#ifndef AFTtv_LN_updates_H
#define AFTtv_LN_updates_H

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

void AFTtv_LN_update_btv(const arma::mat &Ymat,
                         arma::vec &VyL,
                         arma::vec &VyU,
                         arma::vec &Vc0,
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
                         arma::vec &accept_btv);

void AFTtv_LN_update_beta(arma::vec &VyL,
                          arma::vec &VyU,
                          arma::vec &Vc0,
                          const arma::vec &yUInf,
                          const arma::vec &yLUeq,
                          const arma::vec &c0Inf,
                          const arma::mat &Xmat,
                          arma::vec &beta,
                          double &mu,
                          double &sigSq,
                          double &beta_prop_var,
                          arma::vec &accept_beta);

void AFTtv_LN_update_mu(const arma::vec &logVyL,
                        const arma::vec &logVyU,
                        const arma::vec &logVc0,
                        const arma::vec &yUInf,
                        const arma::vec &yLUeq,
                        const arma::vec &c0Inf,
                        double &mu,
                        double &sigSq,
                        double &mu_prop_var,
                        int &accept_mu);

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
                           int &accept_sigSq);

#endif

