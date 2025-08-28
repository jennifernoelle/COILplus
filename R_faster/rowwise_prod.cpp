#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rowwise_prod(NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  NumericVector out(nrow, 1.0);
  
  for (int i = 0; i < nrow; ++i) {
    double acc = 1.0;
    for (int j = 0; j < ncol; ++j) {
      acc *= mat(i, j);
    }
    out[i] = acc;
  }
  
  return out;
}
