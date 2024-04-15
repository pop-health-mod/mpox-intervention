// author: Jesse Knight
#include <Rcpp.h>
#include <boost/multi_array.hpp>
#define TOL 1e-12

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix ipf_cpp(NumericMatrix M, NumericVector m1, NumericVector m2){
  NumericVector r1, r2;
  int i, k;
  for (k = 0; k < 100; k++){
    r1 = m1 / rowSums(M);
    for (i = 0; i < M.nrow(); i++ ){ M(i,_) = M(i,_) * r1(i); }
    r2 = m2 / colSums(M);
    for (i = 0; i < M.ncol(); i++ ){ M(_,i) = M(_,i) * r2(i); }
    if (is_true(all( abs(r1-1) < TOL & abs(r2-1) < TOL))){ break; }
  }
  return M;
}
