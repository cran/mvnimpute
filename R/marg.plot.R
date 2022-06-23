#' Marginal density plots function
#'
#' Draws the marginal density plots for all variables.
#'
#' @param data.mat data matrix including all the variables.
#' @param title title of the generated plots.
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
#' ### marginal density plots of variables in the dataset from the last iteration
#' marg.plot(sim.res$imputed.data[[num.iter]], title = colnames(miss.dat))
#'
#' @return Marginal density plot for each variable in the dataset.
#'
#' @export
marg.plot <- function(data.mat, title = NULL) {

  # draw marginal plots to visual distribution of each marginal
  if (!is.null(ncol(data.mat))) {

    for (i in 1:ncol(data.mat)) {

      plot(density(na.omit(data.mat[, i])), main = title[i])

    }

  } else {

    plot(density(na.omit(data.mat)), main = title)

  }

}
