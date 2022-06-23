#' Averaged simulated values plot function
#'
#' Calculates and the average simulated values of all parameters and generates plots.
#'
#' @param data.mat data matrix including the simulated values for plot.
#' @param start the number of cycle to start.
#' @param end the number of cycle to end.
#' @param x.lab label of the x axis in the generated plot, default is set to "Iteration number".
#' @param y.lab label of the y axis in the generated plot, default is set to "Average of simulated values".
#' @param title title of the generated plot.
#' @param details logical variable to specify whether the average simulated values are returned, default is set to FALSE.
#'
#' @details This function calculates the average simulated values across simulations.
#' \code{iter} can be any number of iterations you want to draw, the corresponding number of rows
#' of the data should be \code{iter} + 1.
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
#' ### plot of average values of simulated mean and variance
#' avg.plot(sim.res$simulated.mu, 201, 501, title = "Simulated value: mean")
#' avg.plot(sim.res$simulated.sig, 201, 501, title = "Simulated value: variance")
#'
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
