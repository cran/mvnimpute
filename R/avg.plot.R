#' Averaged simulated values plot function
#'
#' Calculates the average simulated values of all parameters and generates plots.
#'
#' @param data.mat data matrix including the simulated values for plot.
#' @param start the number of cycle to start.
#' @param end the number of cycle to end.
#' @param x.lab label of the x axis in the generated plot, default is set to "Iteration number".
#' @param y.lab label of the y axis in the generated plot, default is set to "Average of simulated values".
#' @param title title of each generated plot.
#' @param details logical variable to specify whether the average simulated values are returned, default is set to FALSE.
#'
#' @details This function calculates the average simulated values across simulations.
#' \code{iter} can be any number of iterations you want to draw, the corresponding number of rows
#' of the data should be \code{iter} + 1.
#'
#' @examples
#'
#' ### generate some normal data
#' dat <- MASS::mvrnorm(n = 1000, mu = c(1, 2, 3, 4), Sigma = diag(4))
#'
#' ### set column names
#' colnames(dat) <- paste0("Var ", 1:ncol(dat))
#'
#' ### average values plot: take sample from 500 to 1000 rows
#' avg.plot(data.mat = dat[500:1000, ], start = 500, end = 1000, title = "Random Variables")
#'
#' @return The plot of averaged values across iterations. If \code{details} = TRUE,
#' a matrix containing the averaged values of all the variables across iterations will be returned.
#'
#' @export
avg.plot <- function(data.mat,  ### matrix that includes values for plot
                     start,     ### the index of the first iteration for drawing
                     end,       ### the index of the last iteration for drawing
                     x.lab = "Iteration number",              ### label of x axis
                     y.lab = "Average of simulated values",   ### label of y axis
                     title = NULL,
                     details = FALSE) {  ### print the averages of simulated values in the console

  # check column names
  if (is.null(colnames(data.mat))){stop("Variable names have to be specified")}
  # compute the cumulative averages of the parameter values to the current simulation
  iter <- end - start
  avg.param <- matrix(nrow = iter  + 1, ncol = ncol(data.mat))
  avg.param[1, ] <- data.mat[1, ]

  for (i in 2:nrow(avg.param)) {

    avg.param[i, ] <- apply(data.mat[1:i, ], 2, mean)

  }

  # transform data from wide to long

  wide.dat <- as.data.frame(cbind(avg.param, start:end))
  colnames(wide.dat) <- c(colnames(data.mat), "iter")

  # make sure the id column should be factor variable
  # wide.mean$iter <- as.factor(wide.mean$iter)

  long.dat <- melt(wide.dat, id.vars = "iter",
                   measure.vars = colnames(wide.dat)[colnames(wide.dat) != "iter"],
                   variable.name = "var",
                   value.name = "value")

  print(ggplot(data = long.dat,
         aes(x = as.numeric(iter),
             y = .data$value,
             group = var,
             color = var)) +
    geom_line() +
    xlab(x.lab) +
    ylab(y.lab) +
    ggtitle(title))

  if (details) return(avg.param)

}
