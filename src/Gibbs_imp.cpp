// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
using namespace Rcpp;

/***********************************/
// Function for performing Gibbs sampler for the imputation step for multiple imputation
// [[Rcpp::export]]
NumericMatrix Gibbs_imp(const NumericMatrix data, const List data_indx, const NumericVector mu_vec, const NumericMatrix cond_param) {

  const int n = data.nrow(), p = data.ncol();
  // bounds information
  const NumericMatrix ll = data_indx[0];
  const NumericMatrix ul = data_indx[1];

  // preventing the from modifying the data matrix
  NumericMatrix iter_data(data.nrow(), data.ncol());
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      iter_data(i, j) = data(i, j);
    }
  }

  for (int i = 0; i < n; i++) {

    for (int j = 0; j < p; j++) {
      // covariates vector
      NumericVector cov_ij(p - 1);
      for (int k = 0, count = 0; k < p; k++) {
        if (k != j) {
          cov_ij(count) = iter_data(i, k) - mu_vec(k);
          count++;
        }
      }

      // conditional parameters
      NumericVector cond_j = cond_param(j, _);
      double mean_i = 0.0;
      // calculate mean
      for (int k = 0; k < p - 1; k++) {
        mean_i += cond_j(k + 1) * cov_ij(k);
      }

      // (1) missing values
      // if (ll(i, j) != ul(i, j) && (abs(ll(i, j)) == abs(ul(i, j))))
      // change finite censoring limits to NA values -- 6.9.2022
      if (LogicalVector::is_na(ll(i, j)) && LogicalVector::is_na(ul(i, j)))
        iter_data(i, j) = R::rnorm(mu_vec(j) + mean_i, sqrt(cond_j(p)));
      // (2) censored values
      // else if (ll(i, j) != ul(i, j) && (abs(ll(i, j)) != abs(ul(i, j))))
      // change finite censoring limits to NA values -- 6.9.2022
      else if ((!LogicalVector::is_na(ll(i, j)) || !LogicalVector::is_na(ul(i, j))) &&
        (ll(i, j) != ul(i, j))) {
        double lower_val = ll(i, j);
        double upper_val = ul(i, j);
        if (LogicalVector::is_na(ll(i, j)))
          lower_val = -1000000;  // pick very small number for the lower limit
        if (LogicalVector::is_na(ul(i, j)))
          upper_val = 1000000;   // pick very large number for the upper limit
        iter_data(i, j) = r_truncnorm(mu_vec(j) + mean_i, sqrt(cond_j(p)), lower_val, upper_val);
      }
    }
  }

  return iter_data;
}

/***********************************/
