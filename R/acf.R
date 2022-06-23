#' Autocorrelation function
#'
#' Calculates the autocorrelation function and draws the plots.
#'
#' @param data.mat matrix including the variables of which autocorrelations are calculated.
#' @param lag lag at which the autocorrelation is calculated, default is set as 50.
#' @param plot  logical variable to specify whether the plot is generated, default is set to TRUE.
#' @param title title of the generated autocorrelation plots.
#' @param details boolean variable to specify whether the autocorrelation values are returned, default is set to FALSE.
#'
#' @details This function calculates the autocorrelations of all the variables on a column by column base.
#' The default value of \code{lag} is set as 50, the maximum number of lag should not exceed the number of rows of the dataset,
#' which reflects the corresponding number of iteration of running the multiple imputation.
#'
#' @examples
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
#' ### ACF of simulated mean and variance
#' acf.calc(sim.res$simulated.mu, title = "ACF: mean")
#' acf.calc(sim.res$simulated.sig, title = "ACF: variance")
#'
#' @return If \code{details} = TRUE, a matrix containing the calculated autocorrelations of all the variables in the dataset will be returned.
#' If \code{plot} = TRUE, the autocorrelation plots of all the variables will be drawn.
#'
#' @export
acf.calc <- function(data.mat, lag = 50, plot = TRUE, title = NULL, details = FALSE) {

  auto.cor <- matrix(NA, nrow = lag + 1, ncol = ncol(data.mat))

  for (i in 1:ncol(data.mat)) {
    # the first row of data matrix is the initial values
    for (j in 0:lag) auto.cor[j + 1, i] <- round(cor(data.mat[2:(nrow(data.mat) - j), i],
                                                     data.mat[(j + 2):nrow(data.mat), i]), 3)
  }
  rownames(auto.cor) <- paste("lag ", 0:lag, sep = "")
  colnames(auto.cor) <- colnames(data.mat)

  if (plot) {
    for (i in 1:ncol(auto.cor)) {
      plot(auto.cor[, i], type = "h", ylim = c(-0.5, 1.0),
           xlab = "lag", ylab = "acf",
           main = title[i])
      abline(h = 0)
    }
  }

  if (details) return(auto.cor)
}
