// [[Rcpp::depends(RcppDist)]]
#include <RcppDist.h>
using namespace Rcpp;

// [[Rcpp::export]]
List param_calc(const List data) {

  NumericMatrix lval = data[0];
  NumericMatrix rval = data[1];

  const int n = lval.nrow(), p = lval.ncol();

  // CC parameters
  NumericMatrix CC_ind(n, p);
  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      // if (lval(i, j) == rval(i, j))
      // change finite censoring limits to NA values -- 6.9.2022
      if ((!LogicalVector::is_na(lval(i, j)) && !LogicalVector::is_na(rval(i, j))) &&
          lval(i, j) == rval(i, j))
        CC_ind(i, j) = 1;
      else
        CC_ind(i, j) = 0;
    }
  }

  // count CC case
  int CC_num = 0;
  NumericVector CC_row = rowSums(CC_ind);
  for (int i = 0; i < n; i++) {
    if (CC_row(i) == p)
      CC_num++;
  }

  NumericMatrix CC_dat(CC_num, p);
  int ct = 0;
  for (int i = 0; i < n; i++) {
    if (CC_row(i) == p) {
      CC_dat(ct, _) = lval(i, _);
      ct++;
    }
  }

  NumericVector CC_mean = colMeans(CC_dat);
  NumericVector CC_var(p);

  for (int i = 0; i < p; i++) {
    CC_var(i) = var(CC_dat(_, i));
  }

  // AC parameters

  NumericVector not_miss(p);  // count number of observed values
  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      // if (lval(i, j) == rval(i, j))
      // change finite censoring limits to NA values -- 6.9.2022
      if ((!LogicalVector::is_na(lval(i, j)) && !LogicalVector::is_na(rval(i, j))) &&
          lval(i, j) == rval(i, j))
        not_miss(j)++;
    }
  }

  NumericVector AC_mean(p);
  NumericVector AC_var(p);

  // CC imputation
  for (int j = 0; j < p; j++) {
    int len_var = not_miss(j);
    NumericVector obs_val(len_var);

    // observed data
    for (int i = 0, count = 0; i < n; i++) {
      // if (lval(i, j) == rval(i, j)) {
      // change finite censoring limits to NA values -- 6.9.2022
      if ((!LogicalVector::is_na(lval(i, j)) && !LogicalVector::is_na(rval(i, j))) &&
          lval(i, j) == rval(i, j)) {
        obs_val(count) = lval(i, j);
        count++;
      }
    }

    AC_mean(j) = mean(obs_val);
    AC_var(j) = var(obs_val);
  }

  return List::create(
    Named("CC.mean") = CC_mean,
    Named("CC.var") = CC_var,
    Named("AC.mean") = AC_mean,
    Named("AC.var") = AC_var
  );

}

// [[Rcpp::export]]
NumericMatrix single_imputation(List data) {

  NumericVector AC_mean = param_calc(data)["AC.mean"];
  NumericVector AC_var = param_calc(data)["AC.var"];

  NumericMatrix lval = data[0], rval = data[1];

  // number of observations and variables
  const int n = lval.nrow(), p = lval.ncol();

  // SI data
  NumericMatrix SI_dat(n, p);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      // (1) observed values
      // if (lval(i, j) == rval(i, j))
      // change finite censoring limits to NA values -- 6.9.2022
      if (!(LogicalVector::is_na(lval(i, j)) || LogicalVector::is_na(rval(i, j))) &&
          lval(i, j) == rval(i, j))
        SI_dat(i, j) = lval(i, j);
      // (2) missing and censored values
      else
        SI_dat(i, j) = R::rnorm(AC_mean(j), sqrt(AC_var(j)));
      // (3) censored values
      // else {
      //   SI_dat(i, j) = r_truncnorm(CC_mean(j), sqrt(CC_var(j)), lval(i, j), rval(i, j));
      // }
    }
  }

  // return SI_dat;
  return SI_dat;

}
