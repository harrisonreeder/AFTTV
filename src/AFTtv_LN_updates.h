#ifndef AFTtv_LN_updates_H
#define AFTtv_LN_updates_H

double AFTtv_LN_loglik(const arma::vec &logVyL,
                       const arma::vec &logVyU,
                       const arma::vec &logVc0,
                       const arma::vec &logvyL,
                       const arma::vec &yUInf,
                       const arma::vec &yLUeq,
                       const arma::vec &c0Inf,
                       double &mu,
                       double &sigSq);

void AFTtv_LN_update_mu(const arma::vec &logVyL,
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

void AFTtv_LN_update_sigSq(const arma::vec &logVyL,
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

void AFTtv_LN_update_beta(arma::vec &VyL,
                          arma::vec &VyU,
                          arma::vec &Vc0,
                          arma::vec &vyL,
                          const arma::vec &yUInf,
                          const arma::vec &yLUeq,
                          const arma::vec &c0Inf,
                          const arma::mat &Xmat,
                          arma::vec &beta,
                          double &mu,
                          double &sigSq,
                          double &beta_prop_var,
                          arma::vec &accept_beta,
                          double &curr_loglik);


void AFTtv_LN_update_hyper_icar(const arma::vec &beta_tv,
                                double &meanbtv,
                                double &varbtv,
                                double &a_btv,
                                double &b_btv,
                                const arma::mat &cholinvSigma_btv);

void AFTtv_LN_update_btv_icar(const arma::mat &Ymat,
                              arma::vec &logVyL,
                              arma::vec &logVyU,
                              arma::vec &logVc0,
                              arma::vec &logvyL,
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
                              double &meanbtv,
                              double &varbtv,
                              const arma::mat &cholinvSigma_btv,
                              arma::vec &accept_btv,
                              double &curr_loglik);

void AFTtv_LN_update_btv(const arma::mat &Ymat,
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

