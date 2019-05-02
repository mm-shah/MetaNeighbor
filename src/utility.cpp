#include <Rcpp.h>
using namespace Rcpp;

template<int RTYPE>
NumericMatrix normalize_cols_cpp_imp(Matrix<RTYPE> M) {
  NumericMatrix result(M.nrow(), M.ncol());
  for (int j = 0; j < M.ncol(); j++) {
    double mean = 0;
    for (int i = 0; i < M.nrow(); i++) { mean += M(i,j); }
    mean /= M.nrow();
    for (int i = 0; i < M.nrow(); i++) { result(i,j) = M(i,j) - mean; }
    double norm = 0;
    for (int i = 0; i < M.nrow(); i++) { norm += result(i,j) * result(i,j); }
    norm = 1 / sqrt(norm);
    for (int i = 0; i < M.nrow(); i++) { result(i,j) *= norm; }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix normalize_cols_cpp(SEXP M) {
  switch (TYPEOF(M)) {
    case INTSXP: return normalize_cols_cpp_imp<INTSXP>(M);
    case REALSXP: return normalize_cols_cpp_imp<REALSXP>(M);
  }
  return 0;
}

// [[Rcpp::export]]
NumericVector count_to_rank(NumericVector x, int total_count) {
  NumericVector result(x.length());
  int cumulated_total = 0;
  for (int i = 0; i < x.length(); i++) {
    result[i] = (cumulated_total + (x(i)+1)*0.5) / total_count;
    cumulated_total += x[i];
  }
  return result;
}

// [[Rcpp::export]]
void bin_to_rank(NumericVector bins, NumericVector rank_per_bin) {
  for (int i = 0; i < bins.length(); i++) { bins[i] = rank_per_bin[bins[i]-1]; }
}
