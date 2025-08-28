#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector row_logprod_mask_idx(const IntegerMatrix& curr_inter,   // n x m (0/1)
                                   const IntegerMatrix& focus_st,     // n x m (0/1)
                                   const IntegerVector& occ_idx,      // indices of j with occur_others==1
                                   const NumericMatrix& log1m_pipj)   // n x m
{
  const int n = curr_inter.nrow();
  NumericVector out(n);
  const int k = occ_idx.size();
  
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int t = 0; t < k; ++t) {
      const int j = occ_idx[t] - 1;  // R->C++ index
      if (curr_inter(i,j) && focus_st(i,j)) {
        acc += log1m_pipj(i,j);
      }
    }
    out[i] = acc;  // this is log(prod)
  }
  return out;
}
