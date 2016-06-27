


#' Calculate 3D covariance matrix using unscaled congruence coefficient.  Skips any missing values in computation of covariance matrix 
#'
#' Calculate 3D covariance matrix using unscaled congruence coefficient.  
#'   Skips any missing values in computation of covariance matrix
#' 
#'
#'@param M An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@export
#'
#'@details This function does not guarantee that the returned matrix is  
#'  positive definite. If the covariance matrix is not positive definite 
#'  a warning is given and the matrix can be bent to create the closest
#'  positive definite matrix with \code{as.matrix(Matrix::nearPD(mat)$mat)}.
#'
#'@return N x N covariance matrix
#'@examples
#' M <- array(rnorm(4 * 2 * 3), dim = c(2, 3, 4)) 
#' M.cvm <- dotcvm(M)

dotcvm <- function(M){
 # Calculate covariance between each pairs of columns.
  N <- matrix(NA, nrow = dim(M)[1], ncol = dim(M)[1])
  for(i in 1:dim(M)[1]){
    for(j in i:dim(M)[1]){
      N[i, j] <- dotcvmentry(M, i, j)
    }
  }
  
  N[lower.tri(N)] <- t(N)[lower.tri(N)]
  
  e <- min(eigen(N)$values)
  if(e < 0) warning(paste('CVM has negative eigenvalue', e))
  return(N)
}


# Calculate 3D covariance between two landmarks across specimens 
# 
#
#@param M An N x 3 x M array. M = no of specimens, N = no of landmarks.
#@param col1 Integer of first column 
#@param col2 Integer of second column 
#
#@return Covariance value
#@name dotcvmentry



# Check columns have enough data and then calculate covariance between columns
dotcvmentry <- function(M, col1, col2){
  n <- 0
  s1 <- c(0, 0, 0)
  s2 <- c(0, 0, 0)

  # For each specimen
  for(i in 1:dim(M)[3]){
    if(!anyNA(M[c(col1, col2), , i])){
      n <- n + 1
      s1 <- s1 + M[col1, , i]
      s2 <- s2 + M[col2, , i]
    }
  }

  if(n <= 1) stop(paste("There is too much missing data to covary columns", col1, "and", col2))

  s1 <- s1 / n
  s2 <- s2 / n

  p <- 0

  # for each specimen
  for(i in 1:dim(M)[3]){
    if(!anyNA(M[c(col1, col2), , i])){
      p <- p + crossprod((M[col1, , i] - s1), (M[col2, , i] - s2))
    }
  }

  return(p / (n - 1))
}







#' Calculate 3D correlation matrix using the congruence coefficient.  Skips any missing values in computation of correlation matrix 
#' 
#'  Calculate 3D correlation matrix using the congruence coefficient.  
#'    Skips any missing values in computation of correlation matrix
#'    Gives an N x N correlation matrix.
#'
#'@param M An N x D x M array. M = no of specimens, N = no of landmarks, D = 3 dimensions
#'@export
#'
#'@return Correlation matrix
#'@examples
#' M <- array(rnorm(4 * 2 * 3), dim = c(2, 3, 4)) 
#' M.corr <- dotcorr(M)
#'


dotcorr <- function(M){
  # Calculate covariance between each pairs of columns.
  N <- matrix(NA, nrow = dim(M)[1], ncol = dim(M)[1])
  for(i in 1:dim(M)[1]){
    for(j in i:dim(M)[1]){
      N[i, j] <- dotcorrentry(M, i, j)
    }
  }
  
  N[lower.tri(N)] <- t(N)[lower.tri(N)]
  
  e <- min(eigen(N)$values)
  if(e < 0) warning(paste('CVM has negative eigenvalue', e))
  return(N)
}









# Calculate 3D correlation between two landmarks across specimens 
# 
# Calculates the congruence coefficient for 3 dimensional landmarks
#
#@param M An N x 3 x M array. M = no of specimens, N = no of landmarks.
#@param col1 Integer of first column 
#@param col2 Integer of second column 
#
#@return Length 1 numeric giving the congruence correlation



# Check columns have enough data and then calculate correlation between columns
dotcorrentry <- function(M, col1, col2){
  n <- 0
  s1 <- c(0, 0, 0)
  s2 <- c(0, 0, 0)

  # For each specimen
  for(i in 1:dim(M)[3]){
    if(!anyNA(M[c(col1, col2), , i])){
      n <- n + 1
      s1 <- s1 + M[col1, ,i]
      s2 <- s2 + M[col2, , i]
    }
  }

  if(n <= 1) stop(paste("There is too much missing data covary columns", col1, "and", col2))

  s1 <- s1/n
  s2 <- s2/n

  p <- 0
  sumi <- 0
  sumj <- 0

  # for each specimen
  for(i in 1:dim(M)[3]){
    if(!anyNA(M[c(col1, col2), , i])){
      p <- p + crossprod((M[col1, , i] - s1), (M[col2, , i] - s2))
	    sumi <- sumi + crossprod((M[col1, , i] - s1), (M[col1, , i] - s1))
	    sumj <- sumj + crossprod((M[col2, , i] - s2), (M[col2, , i] - s2))
    }
  }

  return(p / sqrt(sumi * sumj))
}









#' Calculate covariance matrix between individual dimensions within landmarks
#'
#' Calculate covariance matrix between individual dimensions within landmarks. Skips any missing values
#'   in computation of covariance matrix.
#' 
#'
#'@param M An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@export
#'
#'@details This function does not guarantee that the returned matrix is  
#'  positive definite. If the covariance matrix is not positive definite 
#'  a warning is given and the matrix can be bent to create the closest
#'  positive definite matrix with \code{as.matrix(Matrix::nearPD(mat)$mat)}.
#'
#'@return 3N x 3N covariance matrix
#'@examples
#' M <- array(rnorm(4 * 2 * 3), dim = c(2, 3, 4)) 
#' M.cvm <- dotcovar(M)

dotcovar <- function(M){
  # Calculate covariance between each pairs of columns.

  Mflat <- matrix(NA, nrow = dim(M)[3], ncol = 3 * dim(M)[1])

  Mflat[, seq(1, 3 * dim(M)[1], by = 3)] <- t(M[, 1, ])
  Mflat[, seq(2, 3 * dim(M)[1], by = 3)] <- t(M[, 2, ])
  Mflat[, seq(3, 3 * dim(M)[1], by = 3)] <- t(M[, 3, ])

  N <- matrix(NA, nrow = 3 * dim(M)[1], ncol = 3 * dim(M)[1])

  
  for(i in 1:NCOL(Mflat)){
    for(j in i:NCOL(Mflat)){
      N[i, j] <- cov(Mflat[, i], Mflat[, j], use = 'complete.obs')
    }
  }
  
  N[lower.tri(N)] <- t(N)[lower.tri(N)]
    
  
  e <- min(eigen(N)$values)
  if(e < 0) warning(paste('CVM has negative eigenvalue', e))
  return(N)
}







