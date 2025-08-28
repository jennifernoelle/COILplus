// fast_loglik.cpp
#include <Rcpp.h>
using namespace Rcpp;

// Branchless, numerically stable: y*c - softplus(c)
// c = logit(pi) + log_prod1, softplus(c) = log1p(exp(-|c|)) + max(c, 0)

// [[Rcpp::export]]
NumericVector fast_loglik_log_cpp(NumericVector pi,
                                  NumericVector log_prod1,
                                  NumericVector y) {
  int n = pi.size();
  if (log_prod1.size() != n || y.size() != n)
    stop("length mismatch");
  NumericVector out(n);
  const double eps = 1e-12;
  
  for (int i = 0; i < n; ++i) {
    double p = pi[i];
    if (p < eps) p = eps; else if (p > 1.0 - eps) p = 1.0 - eps;
    double c = std::log(p) - std::log1p(-p) + log_prod1[i];   // logit(p) + log_prod1
    double abs_c = std::fabs(c);
    double sp = std::log1p(std::exp(-abs_c)) + (c > 0.0 ? c : 0.0); // softplus
    double yi = (y[i] != 0.0) ? 1.0 : 0.0;
    out[i] = yi * c - sp;
  }
  return out;
}

