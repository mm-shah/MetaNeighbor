// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// normalize_cols_cpp
NumericMatrix normalize_cols_cpp(SEXP M);
RcppExport SEXP _MetaNeighbor_normalize_cols_cpp(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(normalize_cols_cpp(M));
    return rcpp_result_gen;
END_RCPP
}
// count_to_rank
NumericVector count_to_rank(NumericVector x, int total_count);
RcppExport SEXP _MetaNeighbor_count_to_rank(SEXP xSEXP, SEXP total_countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type total_count(total_countSEXP);
    rcpp_result_gen = Rcpp::wrap(count_to_rank(x, total_count));
    return rcpp_result_gen;
END_RCPP
}
// bin_to_rank
void bin_to_rank(NumericVector bins, NumericVector rank_per_bin);
RcppExport SEXP _MetaNeighbor_bin_to_rank(SEXP binsSEXP, SEXP rank_per_binSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rank_per_bin(rank_per_binSEXP);
    bin_to_rank(bins, rank_per_bin);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MetaNeighbor_normalize_cols_cpp", (DL_FUNC) &_MetaNeighbor_normalize_cols_cpp, 1},
    {"_MetaNeighbor_count_to_rank", (DL_FUNC) &_MetaNeighbor_count_to_rank, 2},
    {"_MetaNeighbor_bin_to_rank", (DL_FUNC) &_MetaNeighbor_bin_to_rank, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MetaNeighbor(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
