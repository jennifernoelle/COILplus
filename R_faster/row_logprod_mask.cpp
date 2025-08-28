#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector row_logprod_mask(const IntegerMatrix& curr_inter,     // num_obs x num_oth (0/1)
                               const IntegerMatrix& focus_st,       // num_obs x num_oth (0/1)
                               const IntegerVector& occ_other_st,   // length num_oth (0/1)
                               const NumericMatrix& log1m_pipj)     // num_obs x num_oth
{
  const int n = curr_inter.nrow();
  const int m = curr_inter.ncol();
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < m; ++j) {
      if (curr_inter(i,j) && focus_st(i,j) && occ_other_st[j]) {
        acc += log1m_pipj(i,j);     // add log(1 - p_i p_j)
      }
    }
    out[i] = acc;                    // row log-product; exp() later if needed
  }
  return out;
}
