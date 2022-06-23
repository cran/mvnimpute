########################################
### new data generation function 7-1 ###
########################################


#' New data generation function
#'
#' Simulates multivariate normal data with missing and censored values. In this function, missing values will be
#' generated first in the multivariate data, then censored values will be generated for the non-missing data.
#'
#' @param num_ind number of subjects.
#' @param mean_vec mean vectors.
#' @param cov_mat covariance matrix.
#' @param miss_var variables that have missing values.
#' @param miss_mech missing mechanism. "MCAR" or "MAR". Default "MCAR".
#' @param miss_prob missing data probability when missing data is MCAR.
#' @param censor_var variables that have censored values.
#' @param censor_type type of censoring. "interval", "right" or "left. Default "interval".
#' @param censor_param rate parameter of the exponential distribution that the censoring times come from.
#'
#' @examples
#' ### generate a multivariate normal dataset of 2000 sample size
#' ### using the default arguments
#' data.generation()
#'
#' @return A list containing the fully observed data, the observed data,
#' the bounds information of the observed data and the data type indicator matrix.
#'
#' @export
data.generation <- function(
  num_ind = 2000,      # number of subjects
  mean_vec = rnorm(5),    # specify mean vectors
  cov_mat = diag(5),      # specify covariance matrix
  miss_var = c(2, 3),     # index of variables that have missing values
  miss_mech = "MCAR", # missing data mechanism
  miss_prob = c(0.2, 0.4),
  censor_var = 4,   # index of variables that have censored values
  censor_type = "interval",
  censor_param = 0.1 # parameter of the distribution that the censored values come from
) {

  # generate the complete multivariate normal data
  full.data <- MASS::mvrnorm(num_ind, mean_vec, cov_mat)
  incomplete <- full.data
  ## number of variables
  p <- length(mean_vec)

  # a list of two containing the bounds information for missing/censored data
  data.indx <- list()
  data.indx[[1]] <- data.indx[[2]] <- incomplete
  indicator <- matrix(1, nrow = num_ind, ncol = p)

  # variables that are fully observed
  ## index of each variable
  var_indx <- 1:p
  ## variables that are fully observed
  obs_var <- var_indx[!var_indx %in% miss_var & !var_indx %in% censor_var]

  ###############################################################
  ## UPDATE 6.9.2022: change finite censoring limits to NA values
  ###############################################################
  # 1. handle missing values
  if (miss_mech == "MCAR") {
    for (i in 1:length(miss_var)) {
      ## the index of missing variables and missing data probabilities
      miss_v <- miss_var[i]
      miss_p <- miss_prob[i]
      for (j in 1:num_ind) {
        incomplete[j, miss_v] <- ifelse(runif(1) <= miss_p, NA, incomplete[j, miss_v])
        ### data index matrix
        data.indx[[1]][j, miss_v] <- ifelse(is.na(incomplete[j, miss_v]), NA, incomplete[j, miss_v])
        data.indx[[2]][j, miss_v] <- ifelse(is.na(incomplete[j, miss_v]), NA, incomplete[j, miss_v])
        ### data type indicator matrix
        if (is.na(incomplete[j, miss_v])) {indicator[j, miss_v] <- 0}
      }
    }
  }
  else if (miss_mech == "MAR") {
    ## missing values depend on the first fully observed variable
    obs_i <- obs_var[1]
    miss_mar <- exp(0.1 + 0.5 * full.data[, obs_i]) / (1 + exp(0.1 + 0.5 * full.data[, obs_i]))
    for (i in 1:length(miss_var)) {
      miss_v <- miss_var[i]

      for (j in 1:num_ind) {
        incomplete[j, miss_v] <- ifelse(runif(1) >= miss_mar[j], NA, incomplete[j, miss_v])
        ### data index matrix
        data.indx[[1]][j, miss_v] <- ifelse(is.na(incomplete[j, miss_v]), NA, incomplete[j, miss_v])
        data.indx[[2]][j, miss_v] <- ifelse(is.na(incomplete[j, miss_v]), NA, incomplete[j, miss_v])
        ### data type indicator matrix
        if (is.na(incomplete[j, miss_v])) {indicator[j, miss_v] <- 0}
      }
    }
  }

  # 2. handle censored values
  if (censor_type == "interval") {

    for (i in 1:length(censor_var)) {
      ## the index of censored variables
      censor_v <- censor_var[i]
      for (j in 1:num_ind) {
        # generate the censoring times
        censor_time <- rexp(2, censor_param)
        min.val <- min(censor_time); max.val <- max(censor_time)

        incomplete[j, censor_v] <- ifelse((incomplete[j, censor_v] >= min.val & incomplete[j, censor_v] <= max.val),
                                          NaN, incomplete[j, censor_v])
        ### data index matrix
        data.indx[[1]][j, censor_v] <- ifelse(is.nan(incomplete[j, censor_v]),
                                              min.val, incomplete[j, censor_v])
        data.indx[[2]][j, censor_v] <- ifelse(is.nan(incomplete[j, censor_v]),
                                            max.val, incomplete[j, censor_v])
        ### data type indicator matrix
        if (is.nan(incomplete[j, censor_v])) {indicator[j, censor_v] <- 4}
      }
    }
  }
  else if (censor_type == "right") {

    for (i in 1:length(censor_var)) {
      ## the index of the censored variables
      censor_v <- censor_var[i]
      for (j in 1:num_ind) {
        ## generate the censoring times
        censor_time <- rexp(1, censor_param)

        incomplete[j, censor_v] <- ifelse(incomplete[j, censor_v] >= censor_time,
                                          censor_time, incomplete[j, censor_v])
        ### data index matrix
        data.indx[[1]][j, censor_v] <- ifelse(incomplete[j, censor_v] == censor_time,
                                              censor_time, incomplete[j, censor_v])
        data.indx[[2]][j, censor_v] <- ifelse(incomplete[j, censor_v] == censor_time,
                                            NA, incomplete[j, censor_v])
        ### data type indicator matrix
        if (incomplete[j, censor_v] >= censor_time) {indicator[j, censor_v] <- 2}
      }
    }
  }
  else {

    for (i in 1:length(censor_var)) {
      ## the index of the censored variables
      censor_v <- censor_var[i]
      for (j in 1:num_ind) {
        ## generate the censoring times
        censor_time <- rexp(1, censor_param)

        incomplete[j, censor_v] <- ifelse(incomplete[j, censor_v] <= censor_time,
                                          censor_time, incomplete[j, censor_v])
        ### data index matrix
        data.indx[[1]][j, censor_v] <- ifelse(incomplete[j, censor_v] == censor_time,
                                            NA, incomplete[j, censor_v])
        data.indx[[2]][j, censor_v] <- ifelse(incomplete[j, censor_v] == censor_time,
                                              censor_time, incomplete[j, censor_v])
        ### data type indicator matrix
        if (incomplete[j, censor_v] <= censor_time) {indicator[j, censor_v] <- 3}
      }
    }
  }

  colnames(full.data) <- colnames(incomplete) <- colnames(data.indx[[1]]) <- colnames(data.indx[[2]]) <-
    colnames(indicator) <- paste0("y", 1:p)

  return(list(full.data = full.data,
              observe.data = incomplete,
              bounds = data.indx,
              indicator = indicator))

}
