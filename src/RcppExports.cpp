// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// AFT_LN_mcmc
Rcpp::List AFT_LN_mcmc(const arma::mat& Wmat, const arma::vec& wUInf, const arma::vec& wLUeq, const arma::vec& c0Inf, const arma::mat& Xmat, const arma::vec& hyper_vec, const arma::vec& tuning_vec, const arma::vec& start_vec, int n_burnin, int n_sample, int thin);
RcppExport SEXP _AFTTV_AFT_LN_mcmc(SEXP WmatSEXP, SEXP wUInfSEXP, SEXP wLUeqSEXP, SEXP c0InfSEXP, SEXP XmatSEXP, SEXP hyper_vecSEXP, SEXP tuning_vecSEXP, SEXP start_vecSEXP, SEXP n_burninSEXP, SEXP n_sampleSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Wmat(WmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type wUInf(wUInfSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type wLUeq(wLUeqSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c0Inf(c0InfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Xmat(XmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type hyper_vec(hyper_vecSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tuning_vec(tuning_vecSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type start_vec(start_vecSEXP);
    Rcpp::traits::input_parameter< int >::type n_burnin(n_burninSEXP);
    Rcpp::traits::input_parameter< int >::type n_sample(n_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(AFT_LN_mcmc(Wmat, wUInf, wLUeq, c0Inf, Xmat, hyper_vec, tuning_vec, start_vec, n_burnin, n_sample, thin));
    return rcpp_result_gen;
END_RCPP
}
// AFTtv_LN_mcmc
Rcpp::List AFTtv_LN_mcmc(const arma::mat& Ymat, const arma::vec& yUInf, const arma::vec& yLUeq, const arma::vec& c0Inf, const arma::mat& Xmat, const arma::vec& Xvec_tv, const arma::vec& prior_vec_num, const arma::vec& hyper_vec, const arma::vec& tuning_vec, const arma::vec& start_vec, const arma::vec& knots_init, int n_burnin, int n_sample, int thin);
RcppExport SEXP _AFTTV_AFTtv_LN_mcmc(SEXP YmatSEXP, SEXP yUInfSEXP, SEXP yLUeqSEXP, SEXP c0InfSEXP, SEXP XmatSEXP, SEXP Xvec_tvSEXP, SEXP prior_vec_numSEXP, SEXP hyper_vecSEXP, SEXP tuning_vecSEXP, SEXP start_vecSEXP, SEXP knots_initSEXP, SEXP n_burninSEXP, SEXP n_sampleSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Ymat(YmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type yUInf(yUInfSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type yLUeq(yLUeqSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c0Inf(c0InfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Xmat(XmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xvec_tv(Xvec_tvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior_vec_num(prior_vec_numSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type hyper_vec(hyper_vecSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tuning_vec(tuning_vecSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type start_vec(start_vecSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type knots_init(knots_initSEXP);
    Rcpp::traits::input_parameter< int >::type n_burnin(n_burninSEXP);
    Rcpp::traits::input_parameter< int >::type n_sample(n_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(AFTtv_LN_mcmc(Ymat, yUInf, yLUeq, c0Inf, Xmat, Xvec_tv, prior_vec_num, hyper_vec, tuning_vec, start_vec, knots_init, n_burnin, n_sample, thin));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_arma
double dmvnrm_arma(const arma::vec& x, double& mean, double& var, const arma::mat& cholinvSigma_sub, int logd);
RcppExport SEXP _AFTTV_dmvnrm_arma(SEXP xSEXP, SEXP meanSEXP, SEXP varSEXP, SEXP cholinvSigma_subSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double& >::type var(varSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cholinvSigma_sub(cholinvSigma_subSEXP);
    Rcpp::traits::input_parameter< int >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_arma(x, mean, var, cholinvSigma_sub, logd));
    return rcpp_result_gen;
END_RCPP
}
// Vx_pw
arma::vec Vx_pw(const arma::vec& T, const arma::mat& Xmat, const arma::vec& beta, const arma::vec& Xvec_tv, const arma::vec& beta_tv, const arma::vec& knots, const int log_out);
RcppExport SEXP _AFTTV_Vx_pw(SEXP TSEXP, SEXP XmatSEXP, SEXP betaSEXP, SEXP Xvec_tvSEXP, SEXP beta_tvSEXP, SEXP knotsSEXP, SEXP log_outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Xmat(XmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xvec_tv(Xvec_tvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_tv(beta_tvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< const int >::type log_out(log_outSEXP);
    rcpp_result_gen = Rcpp::wrap(Vx_pw(T, Xmat, beta, Xvec_tv, beta_tv, knots, log_out));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_AFTTV_AFT_LN_mcmc", (DL_FUNC) &_AFTTV_AFT_LN_mcmc, 11},
    {"_AFTTV_AFTtv_LN_mcmc", (DL_FUNC) &_AFTTV_AFTtv_LN_mcmc, 14},
    {"_AFTTV_dmvnrm_arma", (DL_FUNC) &_AFTTV_dmvnrm_arma, 5},
    {"_AFTTV_Vx_pw", (DL_FUNC) &_AFTTV_Vx_pw, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_AFTTV(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
