#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector fast_loglik_eta_cpp(NumericVector eta,
                                  NumericVector log_prod1,
                                  IntegerVector y) {
  int n = eta.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double c = eta[i] + log_prod1[i];
    double sp = std::log1p(std::exp(-std::fabs(c))) + (c > 0.0 ? c : 0.0);
    out[i] = (y[i] ? c : 0.0) - sp;
  }
  return out;
}

