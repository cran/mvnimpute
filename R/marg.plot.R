#' Marginal density plots function
#'
#' Draws marginal density plots for all variables
#'
#' @param data.mat data matrix including all the variables.
#' @param title title of each generated plot.
#'
#' @examples
#' ### generate some data
#' dat <- MASS::mvrnorm(n = 1000, mu = c(1, 2, 3, 4), Sigma = diag(4))
#'
#' ### set column names
#' colnames(dat) <- paste0("Var ", 1:ncol(dat))
#'
#' ### marginal plots
#' marg.plot(data.mat = dat, title = paste0("Var", 1:nrow(dat)))
#' @return Marginal density plot for each variable in the dataset.
#'
#' @export
marg.plot <- function(data.mat, title = NULL) {

  # check column names
  if (is.null(colnames(data.mat))){stop("Variable names have to be specified")}

  # draw marginal plots to visual distribution of each marginal
  if (!is.null(ncol(data.mat))) {

    for (i in 1:ncol(data.mat)) {

      plot(density(na.omit(data.mat[, i])), main = title[i])

    }

  } else {

    plot(density(na.omit(data.mat)), main = title)

  }

}
