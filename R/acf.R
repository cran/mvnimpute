#' Autocorrelation function
#'
#' Calculates the autocorrelation function and draws the plots.
#'
#' @param data.mat matrix including the variables of which autocorrelations are calculated.
#' @param lag lag at which the autocorrelation is calculated, default is set as 50.
#' @param plot  logical variable to specify whether the plot is generated, default is set to TRUE.
#' @param title title of each generated autocorrelation plot.
#' @param details boolean variable to specify whether the autocorrelation values are returned, default is set to FALSE.
#'
#' @details This function calculates the autocorrelations of all the variables on a column by column base.
#' The default value of \code{lag} is set as 50, the maximum number of lag should not exceed the number of rows of the dataset,
#' which reflects the corresponding number of iteration of running the multiple imputation.
#'
#' @examples
#'
#' ### generate some data
#' dat <- MASS::mvrnorm(n = 1000, mu = c(1, 2, 3, 4), Sigma = diag(4))
#'
#' ### ACF plots
#' acf.calc(data.mat = dat, title = paste0("Var ", 1:nrow(dat)))
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
