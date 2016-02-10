
#' Calculate 1D covariance matrix 
#' 
#' Each row of M should correspond to one _specimen _.  Each column of
#' M should correspond to one _landmark _.  In the 1D case, each element
#' of M should be either a number or a string such as "?" which indicates
#' missing data.  In the 3D case, each element of M should be either
#' a 3-component vector or a string such as "?" which indicates missing
#' data.
#'
#'@param M An M x N matrix. M = no of specimens, N = no of landmarks.
#'@export
#'
#'@return 1D covariance matrix





#' Calculate 3D covariance matrix 
#' 
#' Each row of M should correspond to one _specimen _.  Each column of
#' M should correspond to one _landmark _.  In the 1D case, each element
#' of M should be either a number or a string such as "?" which indicates
#' missing data.  In the 3D case, each element of M should be either
#' a 3-component vector or a string such as "?" which indicates missing
#' data.
#'
#'@param M An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@export
#'
#'@return 3D covariance matrix

dotcvm <- function(M){
  # Calculate covariance between each pairs of columns.
  N <- outer(1:dim(M)[2], 1:dim(M)[2], dotcvmentry, M = M)
  e <- min(eigen(N)$values)
  if(e < 0) warning(paste('CVM has negative eigenvalue', e))
  return(N)
}



#' Calculate 3D covariance between two landmarks across specimens 
#' 
#'
#'@param M An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@col1 Integer of first column 
#'@col2 Integer of second column 
#'
#'@return 3D covariance matrix



# Check columns have enough data and then calculate covariance between columns
dotcvmentry <- Vectorize(function(M, col1, col2){
  n <- 0
  s1 <- c(0, 0, 0)
  s2 <- c(0, 0, 0)

  # For each specimen
  for(i in 1:dim(M)[1]){
    if(!anyNA(M[, c(col1, col2), ])){
      n <- n + 1
      s1 <- s1 + M[i, col1, ]
      s2 <- s2 + M[i, col1, ]
    }
  }

  if(n <= 1) stop(paste("There is too much missing data  covary columns", col1, "and", col2))

  s1 <- s1/n
  s2 <- s2/n

  p <- 0

  # for each specimen
  for(i in 1:dim(M)[1]){
    if(!anyNA(M[, c(col1, col2), ])){
      p <- p + crossprod((M[i, col1, ] - s1), (M[i, col2, ] - s2))
    }
  }

  return(p/(n - 1))
}, vectorize.args=list('col1', 'col2'))







cvmentry <- Vectorize(function(M, col1, col2){
  n <- 0
  s1 <- 0
  s2 <- 0

  # For each specimen
  for(i in 1:dim(M)[1]){
    if(!anyNA(M[, c(col1, col2)])){
      n <- n + 1
      s1 <- s1 + M[i, col1]
      s2 <- s2 + M[i, col1]
    }
  }

  if(n <= 1) stop(paste("There is too much missing data to covary columns", col1, "and", col2))

  s1 <- s1/n
  s2 <- s2/n

  p <- 0

  # for each specimen
  for(i in 1:dim(M)[1]){
    if(!anyNA(M[, c(col1, col2)])){
      p <- p + (M[i, col1] - s1) * (M[i, col2] - s2)
    }
  }

  return(p/(n - 1))

}, vectorize.args=list('col1', 'col2'))



cvm <- function(M){
  # Calculate covariance between each pairs of columns.
  N <- outer(1:dim(M)[2], 1:dim(M)[2], cvmentry, M = M)
  e <- min(eigen(N)$values)
  if(e < 0) warning(paste('CVM has negative eigenvalue', e))
  return(N)
}






