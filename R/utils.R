################################################
## SWEEP Operator from joint model parameters ##
################################################
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





####################
## SWEEP Operator ##
####################

swp <- function
(
  dat, # data input
  swp # the index that you want to sweep on
  # NOTE: swp = 0 is to sweep the first row
)
{
  dim.dat <- dim(dat) # dimension of the input dataset
  # calculate the augmented matrix
  p <- dim.dat[2] # number of variables
  n <- dim.dat[1] # number of observations
  aug.mat <- matrix(NA, nrow = p + 1, ncol = p + 1)
  aug.mat[1, 1] <- 1
  for (i in 1:p) {
    aug.mat[1,(i + 1)] <- aug.mat[(i + 1), 1] <- mean(dat[,i])
    # calculate values on the margins, which should be the means of each variable
    for (j in 1:p) {
      aug.mat[(i + 1), (j + 1)] <- aug.mat[(j + 1), (i + 1)] <- sum(dat[, i] * dat[, j]) / n
      # calculate the SS-CP matrix
    }
  }

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

  for(i in 0:swp) {

    h.jj <- -1 / swp.mat[(i + 1), (i + 1)]
    h.ij <- swp.mat[(i + 1),-(i + 1)] / swp.mat[(i + 1), (i + 1)]
    mat.jk <- swp.mat[-(i + 1), -(i + 1)] - swp.mat[(i + 1),-(i + 1)] %*% t(swp.mat[(i + 1),-(i + 1)]) / swp.mat[(i + 1), (i + 1)]

    swp.mat[i + 1, i + 1] <- h.jj # the swept element
    swp.mat[i + 1, -(i + 1)] <- swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
    swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  }
  return(swp.mat)
}


############################
## Reverse SWEEP Operator ##
############################

Rev.swp <- function(
  dat,
  swp
){

  p <- ncol(dat)

  swp.mat <- swp(dat, p - 1)
  # matrix that has been swept to the greatest extent

  swp.dim <- dim(swp.mat)
  rows <- swp.dim[1]
  cols <- swp.dim[2]

  h.jj <- NA # the diagonal element
  h.ij <- numeric(rows - 1) # the margin vectors
  mat.jk <- matrix(NA, nrow = (rows - 1), ncol = (cols - 1)) # the SS-CP matrix

  for(i in swp:(p - 1)) {

    h.jj <- -1 / swp.mat[(i + 1), (i + 1)]
    h.ij <- -swp.mat[(i + 1),-(i + 1)] / swp.mat[(i + 1), (i + 1)]
    mat.jk <- swp.mat[-(i + 1), -(i + 1)] - swp.mat[(i + 1),-(i + 1)] %*% t(swp.mat[(i + 1),-(i + 1)]) / swp.mat[(i + 1), (i + 1)]

    swp.mat[i + 1, i + 1] <- h.jj # the swept element
    swp.mat[i + 1, -(i + 1)] <- swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
    swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  }
  return(swp.mat)

}

########################################################################################################
## function to calculate the augmented variance-covariance matrix that will be used in SWEEP Operator ##
########################################################################################################

cal.aug.mat <- function (
  dat # data input
)
{
  dim.dat <- dim(dat) # dimension of the input dataset
  # calculate the augmented matrix
  p <- dim.dat[2] # number of variables
  n <- dim.dat[1] # number of observations
  aug.mat <- matrix(NA, nrow = p + 1, ncol = p + 1)
  aug.mat[1, 1] <- 1
  for (i in 1:p) {
    aug.mat[1,(i + 1)] <- aug.mat[(i + 1), 1] <- mean(dat[,i])
    # calculate values on the margins, which should be the means of each variable
    for (j in 1:p) {
      aug.mat[(i + 1), (j + 1)] <- aug.mat[(j + 1), (i + 1)] <- sum(dat[, i] * dat[, j]) / n
      # calculate the SS-CP matrix
    }
  }
  return(aug.mat)
}

#############################
## stepwise SWEEP Operator ##
#############################

step.swp <- function(
  mat, # input matrix to be swept on
  swp # index of row and column to sweep
)
{
  swp.mat <- mat
  i <- swp

  # steps of SWEEP Operator

  h.jj <- -1 / swp.mat[(i + 1), (i + 1)]
  h.ij <- swp.mat[(i + 1),-(i + 1)] / swp.mat[(i + 1), (i + 1)]
  mat.jk <- swp.mat[-(i + 1), -(i + 1)] - swp.mat[(i + 1),-(i + 1)] %*% t(swp.mat[(i + 1),-(i + 1)]) / swp.mat[(i + 1), (i + 1)]

  swp.mat[i + 1, i + 1] <- h.jj # the swept element
  swp.mat[i + 1, -(i + 1)] <- swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
  swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  return(swp.mat)
}

#####################################
## stepwise Reverse SWEEP Operator ##
#####################################


step.rev.swp <- function(
  mat, # inoput matrix to be reversely swept on
  rev.swp # index of row and column to reversely swept
)
{
  rev.swp.mat <- mat
  i <- rev.swp

  # steps of Reverse SWEEP Operator

  h.jj <- -1 / rev.swp.mat[(i + 1), (i + 1)]
  h.ij <- -rev.swp.mat[(i + 1),-(i + 1)] / rev.swp.mat[(i + 1), (i + 1)]
  mat.jk <- rev.swp.mat[-(i + 1), -(i + 1)] - rev.swp.mat[(i + 1),-(i + 1)] %*% t(rev.swp.mat[(i + 1),-(i + 1)]) / rev.swp.mat[(i + 1), (i + 1)]

  rev.swp.mat[i + 1, i + 1] <- h.jj # the swept element
  rev.swp.mat[i + 1, -(i + 1)] <- rev.swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
  rev.swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  return(rev.swp.mat)
}

############################
## conditional parameters ##
############################


conditional.parameters <- function(
  dat
) {

  p <- ncol(dat) # number of variables

  swp.para <- matrix(NA, nrow = p, ncol = p + 1)

  for (i in 1:p) {
    swp.dat <- cbind(dat[, -i], dat[, i]) # i is the index of the variable to be swept on
    swp.para[i, ] <- swp(swp.dat, p - 1)[, p + 1]
  }
  return(swp.para)
}

