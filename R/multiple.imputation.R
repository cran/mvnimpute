#' Multiple imputation function
#'
#' Multiply imputes the missing and censored values in multivariate data.
#'
#'
#' @param data a list of data containing the lower and upper bounds information for the missing and censored values.
#' @param prior.params list of prior parameter specifications.
#' @param initial.values list of initial values.
#' @param iter number of rounds for doing multiple imputation.
#' @param verbose boolean variable indicating whether the running status is printed in the console. Default is set to TRUE.
#'
#' @return A list including the simulated mean and variance values of the assumed normal model, the covariance matrix, the imputed data,
#' and the conditional model parameters across different iterations of multiple imputation.
#'
#' @details A multivariate normal model is assumed on the data, the sweep operator is adopted to
#' calculate the parameters of the conditional models. The implemented multiple imputation algorithm
#' is based on the data augmentation algorithm proposed by Tanner and Wong (1987). The Gibbs sampling algorithm
#' is adopted to update the model parameters and draw imputations of the coarse data. Output is a
#' list including the parameters of the normal models and the imputed data across different iterations
#' of multiple imputation.
#'
#' @examples
#'
#' ### data and indicator
#' miss.dat <- simulated.dat[[1]]
#' data.ind <- simulated.dat[[2]]
#'
#' ### number of observations and variables
#' n <- nrow(miss.dat); p <- ncol(miss.dat)
#'
#' #### bound matrices
#' b1 <- b2 <- matrix(nrow = nrow(data.ind), ncol = ncol(data.ind))
#'
#' for (i in 1:nrow(b1)) {
#'   for (j in 1:ncol(b1)) {
#'     b1[i, j] <- ifelse(data.ind[i, j] != 1, NA,
#'                        miss.dat[i, j])
#'     b2[i, j] <- ifelse(data.ind[i, j] == 0, NA, miss.dat[i, j])
#'   }
#' }
#' colnames(b1) <- colnames(b2) <- colnames(miss.dat)
#'
#' #### create a matrix for including the lower and upper bounds
#' bounds <- list()
#' bounds[[1]] <- b1; bounds[[2]] <- b2
#'
#' ### prior specifications
#' prior.param <- list(
#'   mu.0 = rep(0, p),
#'   Lambda.0 = diag(100, p),
#'   kappa.0 = 2,
#'   nu.0 = p * (p + 1) / 2
#' )
#'
#' ### starting values
#' start.vals <- list(
#'   mu = rep(0, p),
#'   sigma = diag(100, p)
#' )
#'
#'
#' ### MI
#' num.iter <- 500
#'
#' begin <- Sys.time()
#' sim.res <- multiple.imputation(
#'   bounds,
#'   prior.param,
#'   start.vals,
#'   num.iter,
#'   FALSE
#' )
#'
#' @references
#'
#' Goodnight, J. H. (1979). A tutorial on the SWEEP operator. \emph{The American Statistician}, \bold{33(3)}, 149-158.
#'
#' Tanner, M., & Wong, W. (1987). The Calculation of Posterior Distributions by Data Augmentation.
#' \emph{Journal of the American Statistical Association}, \bold{82(398)}, 528-540.
#'
#' @export
multiple.imputation <- function(
  data,           # the list that contains the cnesoring bounds: LL and UL
  prior.params,   # prior specifications
  initial.values, # initial values
  iter,           # iterations of Gibbs sampler
  verbose = TRUE  # boolean variable to print out running status
) {

  ### control statements
  #### check for input: should be a list
  if (!is.list(data)) {
    stop("The input argument should be a list of length two including the lower and upper bounds of the data!")
  }
  #### check for dimension
  if (sum(dim(data[[1]]) == dim(data[[2]])) != 2) {
    stop("Dimension of the two matrices should be equal to the dimension of the observed data!")
  }
  ### check for element of list: should be a matrix
  if (!(is.matrix(data[[1]]) & is.matrix(data[[2]]))) {
    stop("Each element of the input data list should be a matrix!")
  }

  ### single imputation to make up incomplete data
  iter.data <- single_imputation(data)

  # number of observations and variables
  n <- nrow(data[[1]]); p <- ncol(data[[1]])

  ###########################################
  ### prior specification

  #### check for input: should be a list
  if (!is.list(prior.params)) {
    stop("The input argument should be a list of length four including four prior specifications for the NIW prior!")
  }

  mu.0 <- prior.params$mu.0
  Lambda.0 <- prior.params$Lambda.0
  kappa.0 <- prior.params$kappa.0
  nu.0 <- prior.params$nu.0

  #### check for prior specifications
  if(length(mu.0) != p) {
    stop("The length of the prior mean vector should equal to the nubmer of variables!")
  }
  if (sum(dim(Lambda.0) == p) != 2) {
    stop("The prior scale matrix should be a square matrix!")
  }
  if (length(kappa.0) != 1) {
    stop("The prior number of measurements should be a scaler value!")
  }
  if (length(nu.0) != 1) {
    stop("The degrees of freedom should be a scaler value!")
  }

  # initial values
  if (!is.list(initial.values)) {
    stop("This input argument should be a list of length two including the initial values for the mean vector and covariance matrix!")
  }

  mu.iter <- mu.ini <- initial.values$mu
  sig.iter <- sig.ini <- initial.values$sigma

  #### check for initial values
  if (length(mu.iter) != p) {
    stop("The length of the mean vector should equal to the nubmer of variables!")
  }
  if (sum(dim(sig.iter) == p) != 2) {
    stop("The covariance matrix should be a square matrix which has a dimension correspond to the number of variables!")
  }

  # vector and list to store results
  impute <- list()
  Mu.iter <- matrix(nrow = iter + 1, ncol = ncol(iter.data))
  Sig.iter <- matrix(nrow = iter + 1, ncol = ncol(iter.data))
  Covmat <- list()

  ### set the first set of values as the staring values
  Mu.iter[1, ] <- mu.ini
  Sig.iter[1, ] <- diag(sig.ini)
  Covmat[[1]] <- sig.ini
  cond. <- list()

  # rename
  colnames(Mu.iter) <- colnames(Sig.iter) <- colnames(data[[1]])

  # posterior parameters that do not depend on the data
  kappa.n <- kappa.0 + n
  nu.n <- nu.0 + n

  ## Gibbs iteration
  for (i in 1:iter) {

    ##############################################
    ## 1. Use updated parameters to update data ##
    ##############################################

    # SWP to calculate conditional parameters
    cond.param <- cond_param(iter.data)

    # rename the cond.param
    colnames(cond.param) <- c(paste0("beta", 0:(p-1)), "sigma^2")
    rownames(cond.param) <- c("Outcome: y", paste0("Outcome: x", 1:(p-1)))

    ##### I-step
    iter.data <- Gibbs_imp(iter.data, data, mu.iter, cond.param)

    ##############################################
    ## 2. Use updated data to update parameters ##
    ##############################################

    ### posterior parameters
    y.bar <- colMeans(iter.data)
    mu.n <- kappa.0 * mu.0 / (kappa.0 + n) + n  * y.bar / (kappa.0 + n)

    ###### P-step
    # update mu vector from normal distribution condition on Sigma
    mu.iter <- mvrnorm(1, mu.n, sig.iter/kappa.n)

    # posterior parameters for covariance matrix
    S <- t(sweep(iter.data, 2, y.bar)) %*%
      sweep(iter.data, 2, y.bar)

    Lambda.n <- Lambda.0 + S + kappa.0 * n * (y.bar - mu.0) %*% t(y.bar - mu.0) / (kappa.0 + n)

    # update Sigma from inverse-Wishart distribution
    sig.iter <- rinvwishart(nu.n, Lambda.n)

    # rename
    rownames(sig.iter) <- colnames(sig.iter) <- colnames(data[[1]])

    impute[[i]] <- iter.data                                      # store the imputed dataset for i-th iteration
    Mu.iter[i + 1, ] <- mu.iter                                   # store the simulated means from Gibbs sampler
    Sig.iter[i + 1, ] <- diag(sig.iter)                           # store the simulated variances from Gibbs sampler
    Covmat[[i + 1]] <- sig.iter                                   # store the simulated covariance matrices from Gibbs sampler
    cond.[[i]] <- cond.param                                      # store the conditional parameters from SWEEP operator
    if (verbose) message(paste(i, "-th iteration!", sep = ""))    # print out the running status
  }

  ### rename rows of the matrices for mean and variance
  rownames(Mu.iter) <- rownames(Sig.iter) <- c("Initial", paste0("Iteration: ", 1:iter))

  return(list(
    simulated.mu = Mu.iter,          # simulated mean vector: a vector
    simulated.sig = Sig.iter,        # simulated variance vector: a vector
    simulated.cov = Covmat,          # simulate covariance matrix: a list
    imputed.data = impute,           # simulated data: a list
    conditional.params = cond.       # conditional parameters: a list
  ))

}
