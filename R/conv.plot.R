#' @import ggplot2
NULL
#' Convergence plot function
#'
#' Draws convergence plot for the simulated parameter values of all variables.
#'
#' @param data.mat data matrix including the simulated values.
#' @param start the number of cycle to start.
#' @param end the number of cycle to end.
#' @param x.lab label of the x axis in the generated plot, default is set to "Iteration number".
#' @param y.lab label of the y axis in the generated plot, default is set to "Simulated values".
#' @param title title of each generated plot.
#'
#' @examples
#'
#' ### generate some data
#' dat <- MASS::mvrnorm(n = 1000, mu = c(1, 2, 3, 4), Sigma = diag(4))
#'
#' ### set column names
#' colnames(dat) <- paste0("Var ", 1:ncol(dat))
#'
#' ### convergence plot: select samples from 500 to 1000 rows
#' conv.plot(data.mat = dat[500:1000, ], start = 500, end = 1000, title = "Random Variables")
#'
#' @details The function generates the trace plot of simulated values across iterations.
#' \code{iter} can be any number of iterations you want to draw, the corresponding number of rows
#' of the data is \code{iter} + 1.
#'
#' @return The plot of simulated values across iterations.
#'
#' @export
conv.plot <- function(data.mat, ### matrix that includes values for plot
                      start,    ### the index of the first iteration for drawing
                      end,      ### the index of the last iteration for drawing
                      x.lab = "Iteration number",  ### label of x axis
                      y.lab = "Simulated values",  ### label of y axis
                      title = NULL) {

  ## ensure start and end to be logical numbers
  if (end <= start) {stop("end should be bigger than start!")}
  if (!is.numeric(end) | !is.numeric(start)) {stop("Please enter two scalar values!")}
  if ((start < 0) | (end < 0)) {
    stop("Please enter two logical values!")
  }
  # check column names
  if (is.null(colnames(data.mat))){stop("Variable names have to be specified")}

  # transform data from wide to long

  wide.dat <- as.data.frame(cbind(data.mat, start:end))
  colnames(wide.dat) <- c(colnames(data.mat), "iter")

  # make sure the id column should be factor variable

  long.dat <- melt(wide.dat, id.vars = "iter",
                    measure.vars = colnames(wide.dat)[colnames(wide.dat) != "iter"],
                    variable.name = "var",
                    value.name = "value")

  ggplot(data = long.dat,
         aes(x = as.numeric(.data$iter),
             y = .data$value,
             group = var,
             color = var)) +
    geom_line() +
    xlab(x.lab) +
    ylab(y.lab) +
    ggtitle(title)

}
