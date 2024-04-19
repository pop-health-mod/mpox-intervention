// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// matrix_mult_cpp
Rcpp::NumericMatrix matrix_mult_cpp(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2);
RcppExport SEXP _MpoxModelPack_matrix_mult_cpp(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_mult_cpp(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// matrix_round_cpp
NumericMatrix matrix_round_cpp(NumericMatrix mat);
RcppExport SEXP _MpoxModelPack_matrix_round_cpp(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_round_cpp(mat));
    return rcpp_result_gen;
END_RCPP
}
// matrix_exp_cpp
Rcpp::NumericMatrix matrix_exp_cpp(Rcpp::NumericMatrix m1);
RcppExport SEXP _MpoxModelPack_matrix_exp_cpp(SEXP m1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type m1(m1SEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_exp_cpp(m1));
    return rcpp_result_gen;
END_RCPP
}
// outer_cpp
Rcpp::NumericMatrix outer_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y);
RcppExport SEXP _MpoxModelPack_outer_cpp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(outer_cpp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// ipf_cpp
Rcpp::NumericMatrix ipf_cpp(Rcpp::NumericMatrix M, Rcpp::NumericVector m1, Rcpp::NumericVector m2, double tol);
RcppExport SEXP _MpoxModelPack_ipf_cpp(SEXP MSEXP, SEXP m1SEXP, SEXP m2SEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(ipf_cpp(M, m1, m2, tol));
    return rcpp_result_gen;
END_RCPP
}
// apply_mix_odds_cpp
Rcpp::NumericMatrix apply_mix_odds_cpp(Rcpp::NumericMatrix M0, Rcpp::NumericMatrix OR, double tol);
RcppExport SEXP _MpoxModelPack_apply_mix_odds_cpp(SEXP M0SEXP, SEXP ORSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type M0(M0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type OR(ORSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_mix_odds_cpp(M0, OR, tol));
    return rcpp_result_gen;
END_RCPP
}
// fn_model_cpp
Rcpp::List fn_model_cpp(double bbeta_city, double omega_city, double RR_H_city, double RR_L_city, double gamma1_city, bool TRACING, bool VACCINATING);
RcppExport SEXP _MpoxModelPack_fn_model_cpp(SEXP bbeta_citySEXP, SEXP omega_citySEXP, SEXP RR_H_citySEXP, SEXP RR_L_citySEXP, SEXP gamma1_citySEXP, SEXP TRACINGSEXP, SEXP VACCINATINGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type bbeta_city(bbeta_citySEXP);
    Rcpp::traits::input_parameter< double >::type omega_city(omega_citySEXP);
    Rcpp::traits::input_parameter< double >::type RR_H_city(RR_H_citySEXP);
    Rcpp::traits::input_parameter< double >::type RR_L_city(RR_L_citySEXP);
    Rcpp::traits::input_parameter< double >::type gamma1_city(gamma1_citySEXP);
    Rcpp::traits::input_parameter< bool >::type TRACING(TRACINGSEXP);
    Rcpp::traits::input_parameter< bool >::type VACCINATING(VACCINATINGSEXP);
    rcpp_result_gen = Rcpp::wrap(fn_model_cpp(bbeta_city, omega_city, RR_H_city, RR_L_city, gamma1_city, TRACING, VACCINATING));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MpoxModelPack_matrix_mult_cpp", (DL_FUNC) &_MpoxModelPack_matrix_mult_cpp, 2},
    {"_MpoxModelPack_matrix_round_cpp", (DL_FUNC) &_MpoxModelPack_matrix_round_cpp, 1},
    {"_MpoxModelPack_matrix_exp_cpp", (DL_FUNC) &_MpoxModelPack_matrix_exp_cpp, 1},
    {"_MpoxModelPack_outer_cpp", (DL_FUNC) &_MpoxModelPack_outer_cpp, 2},
    {"_MpoxModelPack_ipf_cpp", (DL_FUNC) &_MpoxModelPack_ipf_cpp, 4},
    {"_MpoxModelPack_apply_mix_odds_cpp", (DL_FUNC) &_MpoxModelPack_apply_mix_odds_cpp, 3},
    {"_MpoxModelPack_fn_model_cpp", (DL_FUNC) &_MpoxModelPack_fn_model_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_MpoxModelPack(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
