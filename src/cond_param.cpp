// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix SWP(const NumericMatrix data, const int swp_indx) {

  // dimension of the data
  const int n = data.nrow(), p = data.ncol();
  NumericMatrix aug_mat(p + 1, p + 1);

  aug_mat(0, 0) = 1.0;

  for (int i = 0; i < p; i++) {
    // calculate values on the margins, which should be the means of each variable
    aug_mat(0, i + 1) = mean(data(_, i));
    aug_mat(i + 1, 0) = mean(data(_, i));

    for (int j = 0; j < p; j++) {

      // calculate y_ki * y_kj
      double cross_prod = 0.0;
      for (int k = 0; k < n; k++) {
        cross_prod += data(k, i) * data(k, j);
      }

      aug_mat(i + 1, j + 1) = cross_prod / n;
      aug_mat(j + 1, i + 1) = cross_prod / n;

    }
  }

  // dimension of the augmented covariance matrix
  const int aug_row = aug_mat.nrow(), aug_col = aug_mat.ncol();

  // initiate SWEEP operator
  // create a matrix to store the values after SWEEP operator
  // this matrix will iterate in the loop
  double h_jj = 0.0;
  NumericVector h_ij(aug_row - 1);  // marginal vectors
  NumericMatrix mat_jk(aug_row - 1, aug_col - 1); // SS-CP matrix
  // matrix after sweeping should have same dimension as the augmented matrix
  NumericMatrix swp_mat(aug_row, aug_col);
  for (int i = 0; i < aug_row; i++) {
    for (int j = 0; j < aug_col; j++) {
      swp_mat(i, j) = aug_mat(i, j);
    }
  }

  for (int i = 0; i < swp_indx + 1; i++) {
    // index (i + 1, i + 1)
    h_jj = -1.0 / swp_mat(i, i);

    // marginal vectors
    NumericVector num_ij(aug_row - 1);
    for (int k = 0, count = 0; k < aug_row ; k++) {
      if (k != i) {
        num_ij(count) = swp_mat(i, k);
        count++;
      }
    }
    // index (i + 1, -(i + 1))
    for (int count = 0; count < h_ij.length(); count++) {
      h_ij(count) = num_ij(count) / swp_mat(i, i);
    }
    NumericMatrix mat_prod(aug_row - 1, aug_col - 1);
    for (int k = 0; k < aug_row - 1; k++) {
      for (int j = 0; j < aug_col - 1; j++) {
        mat_prod(k, j) = num_ij(k) * (num_ij(j)) / swp_mat(i, i);
      }
    }

    // index (-(i + 1), -(i + 1))
    NumericMatrix mat_jkp(aug_row - 1, aug_col - 1);

    for (int k = 0, row_ct = 0; k < aug_row; k++) {
      if (k != i) {
        for (int l = 0, col_ct = 0; l < aug_col; l++) {
          if (l != i) {
            mat_jkp(row_ct, col_ct) = swp_mat(k, l);
            col_ct++;
          }
        }
        row_ct++;
      }
    }
    for (int k = 0; k < aug_row - 1; k++) {
      for (int j = 0; j < aug_col - 1; j++) {
        mat_jk(k, j) = mat_jkp(k, j) - mat_prod(k, j);
      }
    }
    /*******************************/
    // swept matrix
    swp_mat(i, i) = h_jj;

    // set marginal vectors
    // row vectors
    for (int j = 0, ct = 0; j < aug_col; j++) {
      if (j != i) {
        swp_mat(i, j) = h_ij(ct);
        ct++;
      }
    }
    // column vectors
    for (int j = 0, ct = 0; j < aug_row; j++) {
      if (j != i) {
        swp_mat(j, i) = h_ij(ct);
        ct++;
      }
    }

    // SS-CP matrix
    for (int k = 0, row_ct = 0; k < aug_row; k++) {
      if (k != i) {
        for (int l = 0, col_ct = 0; l < aug_col; l++) {
          if (l != i) {
            swp_mat(k, l) = mat_jk(row_ct, col_ct);
            col_ct++;
          }
        }
        row_ct++;
      }
    }
  }

  return swp_mat;
}

// calculate conditional parameters
// [[Rcpp::export]]
NumericMatrix cond_param(const NumericMatrix data) {

  // number of variables
  const int n = data.nrow();
  const int p = data.ncol();

  NumericMatrix swp_param(p, p + 1);
  NumericMatrix swp_data(n, p);
  NumericMatrix cond_p(p + 1, p + 1);


  for (int i = 0; i < p; i++) {
    swp_data(_, p - 1) = data(_, i);
    for (int j = 0, ct = 0; j < p; j++) {
      if (j != i) {
        swp_data(_, ct) = data(_, j);
        ct++;
      }
    }
    cond_p = SWP(swp_data, p - 1);

    for (int j = 0; j < p + 1; j++) {
      swp_param(i, j) = cond_p(j, p);
    }
  }

  return swp_param;
}
