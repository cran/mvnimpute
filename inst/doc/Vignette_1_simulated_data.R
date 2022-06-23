## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
old <- options(scipen = 999)
set.seed(115)

## ----setup--------------------------------------------------------------------
library(mvnimpute)

## -----------------------------------------------------------------------------
#### data generation
n <- 500; p <- 6
# set.seed(133)
m <- c(1, 1, 2, 4, 3, 5)
v <- clusterGeneration::genPositiveDefMat(p, "eigen")$Sigma

miss.var <- c(1, 2)
censor.var <- c(3, 5)

example.data <- data.generation(num_ind = n,
                                mean_vec = m,
                                cov_mat = v,
                                miss_var = miss.var,
                                miss_mech = "MAR",
                                miss_prob = NULL,
                                censor_var = censor.var,
                                censor_type = "interval",
                                censor_param = 0.5)

names(example.data)

## -----------------------------------------------------------------------------
# data
full <- example.data$full.data
observed <- example.data$observe.data
censoring.bounds <- example.data$bounds
indicator <- example.data$indicator
# summary
tail(observed) # observed data
tail(censoring.bounds[[1]]);tail(censoring.bounds[[2]]) # censoring bounds information
tail(indicator) # indicator matrix

## ---- message = FALSE, fig.align = "center", fig.height = 5, fig.width = 8----
visual.plot(indicator, title = NULL)

## ---- message = FALSE---------------------------------------------------------
### prior specifications
prior.spec <- list(
  mu.0 = rep(0, p),
  Lambda.0 = diag(10, p),
  kappa.0 = 100,
  nu.0 = p * (p + 1) / 2
)

start.vals <- list(
  mu = rep(1, p),
  sigma = diag(p)
)
### MCMC simulation
iter <- 1000
sim.res <- multiple.imputation(
  censoring.bounds,
  prior.spec,
  start.vals,
  iter,
  TRUE
)

## ---- fig.align = "center", fig.height = 5, fig.width = 5---------------------
conv.plot(sim.res$simulated.mu, 0, iter, title = "convergence plot of the mean values")
conv.plot(sim.res$simulated.sig, 0, iter, title = "convergence plot of the variance values")

## -----------------------------------------------------------------------------
title = paste("acf: variable", 1:p)
acf.calc(sim.res$simulated.mu, title = title)
acf.calc(sim.res$simulated.sig, title = title)

## -----------------------------------------------------------------------------
for (i in c(1, 2, 3, 5)) {
  plot(density(full[, i]), main = colnames(full)[i])
  lines(density(observed[!is.na(observed[, i]), i]), col = 6, lty = 2)
  lines(density(sim.res$imputed.data[[iter]][, i]), col = 4, lty = 3)
}

## ---- echo = FALSE------------------------------------------------------------
##############################################
swp.true <- function(
  aug.mat, # augmented covariance matrix
  # NOTE: swp = 0 is to sweep the first row
  outcome, # outcome variable
  swp.indx # index of variables
) {

  aug.mat.dim <- dim(aug.mat)
  rows <- aug.mat.dim[1]
  cols <- aug.mat.dim[2]

  # initiate SWEEP Operator
  swp.mat <- aug.mat
  # create a matrix to store the values after SWEEP Operator
  # this matrix will iterate in the loop

  h.jj <- NA # the diagonal element
  h.ij <- numeric(rows - 1) # the margin vectors
  mat.jk <- matrix(NA, nrow = (rows - 1), ncol = (cols - 1)) # the SS-CP matrix

  swp.ind <- swp.indx[swp.indx != outcome]

  for (i in swp.ind) {

    h.jj <- -1 / swp.mat[(i + 1), (i + 1)]
    h.ij <- swp.mat[(i + 1),-(i + 1)] / swp.mat[(i + 1), (i + 1)]
    mat.jk <- swp.mat[-(i + 1), -(i + 1)] - swp.mat[(i + 1),-(i + 1)] %*% t(swp.mat[(i + 1),-(i + 1)]) / swp.mat[(i + 1), (i + 1)]

    swp.mat[i + 1, i + 1] <- h.jj # the swept element
    swp.mat[i + 1, -(i + 1)] <- swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
    swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  }
  return(swp.mat)
}

# true regression parameters
## construct the augmented covariance matrix
aug.cov <- rbind(c(-1, m), cbind(m, v))
colnames(aug.cov) <- rownames(aug.cov) <- NULL
out.var <- 4
true.beta <- swp.true(aug.cov, out.var, 1:p)[, out.var + 1][-(out.var + 1)]

## -----------------------------------------------------------------------------
full.reg <- coef(lm(y4 ~ y1 + y2 + y3 + y5 + y6, data = data.frame(full)))
cc.reg <- coef(lm(y4 ~ y1 + y2 + y3 + y5 + y6, data = data.frame(observed)))

### imputed data
gibbs <- seq(800, 1000, 50)

#### mvnimpute
mvnimpute.dat <- list()
model.param <- list()
for (i in 1:length(gibbs)) {
  
  mvnimpute.dat[[i]] <- sim.res$imputed.data[[gibbs[i]]]
  colnames(mvnimpute.dat[[i]]) <- paste0("y", 1:p)
  model.param[[i]] <- lm(y4 ~ y1 + y2 + y3 + y5 + y6, data.frame(mvnimpute.dat[[i]]))
  
}

mvnimpute.mod <- summary(mice::pool(model.param))
sim.reg <- mvnimpute.mod[, 2]

reg.compare <- data.frame(
  true = true.beta,
  full = full.reg,
  mvnimpute = sim.reg,
  cc = cc.reg
)
reg.compare

options(old)

