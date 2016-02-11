


#' Calculate 3D covariance matrix 
#' 
#'
#'@param M An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@export
#'
#'@return N x N covariance matrix
#'@name dotcvmd

dotcvm <- function(M){
 # Calculate covariance between each pairs of columns.
  N <- matrix(NA, nrow = dim(M)[2], ncol = dim(M)[2])
  for(i in 1:dim(M)[2]){
    for(j in i:dim(M)[2]){
      N[i, j] <- dotcorrentry(M, i, j)
    }
  }
  
  N[lower.tri(N)] <- t(N)[lower.tri(N)]
  
  e <- min(eigen(N)$values)
  if(e < 0) warning(paste('CVM has negative eigenvalue', e))
  return(N)
}


#' Calculate 3D covariance between two landmarks across specimens 
#' 
#'
#'@param M An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@param col1 Integer of first column 
#'@param col2 Integer of second column 
#'
#'@return Covariance value
#'@name dotcvmentry



# Check columns have enough data and then calculate covariance between columns
dotcvmentry <- function(M, col1, col2){
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

  if(n <= 1) stop(paste("There is too much missing data to covary columns", col1, "and", col2))

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
}






#' Calculate 1D covariance between two landmarks across specimens 
#' 
#'
#'@param M An M x N matrix. M = no of specimens, N = no of landmarks.
#'
#'@return 1D covariance matrix



cvm <- function(M){
  # Calculate covariance between each pairs of columns.
  N <- outer(1:dim(M)[2], 1:dim(M)[2], cvmentry, M = M)
  e <- min(eigen(N)$values)
  if(e < 0) warning(paste('CVM has negative eigenvalue', e))
  return(N)
}


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

  if(n <= 1) stop(paste("There is too much missing data  covary columns", col1, "and", col2))

  s1 <- s1/n
  s2 <- s2/n

  p <- 0

  # for each specimen
  for(i in 1:dim(M)[1]){
    if(!anyNA(M[, c(col1, col2)])){
      p <- p + (M[i, col1] - s1) * (M[i, col2] - s2)
    }
  }

  return(p / (n - 1))

}, vectorize.args=list('col1', 'col2'))




#' Calculate 3D correlation matrix 
#' 
#' Calculates the congruence coefficient for 2 or 3 dimensional landmarks
#'   to give a M x M correlation matrix.
#'
#'@param M An M x N x D array. M = no of specimens, N = no of landmarks, D = 2 or 3 dimensions
#'@export
#'
#'@return Correlation matrix


dotcorr <- function(M){
  # Calculate covariance between each pairs of columns.
  N <- matrix(NA, nrow = dim(M)[2], ncol = dim(M)[2])
  for(i in 1:dim(M)[2]){
    for(j in i:dim(M)[2]){
      N[i, j] <- dotcorrentry(M, i, j)
    }
  }
  
  N[lower.tri(N)] <- t(N)[lower.tri(N)]
  
  e <- min(eigen(N)$values)
  if(e < 0) warning(paste('CVM has negative eigenvalue', e))
  return(N)
}









#' Calculate 3D correlation between two landmarks across specimens 
#' 
#' Calculates the congruence coefficient for 2 or 3 dimensional landmarks
#'
#'@param M An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@param col1 Integer of first column 
#'@param col2 Integer of second column 
#'
#'@return Length 1 numeric giving the congruence correlation



# Check columns have enough data and then calculate correlation between columns
dotcorrentry <- function(M, col1, col2){
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
  sumi <- 0
  sumj <- 0

  # for each specimen
  for(i in 1:dim(M)[1]){
    if(!anyNA(M[, c(col1, col2), ])){
      p <- p + crossprod((M[i, col1, ] - s1), (M[i, col2, ] - s2))
	    sumi <- sumi + crossprod((M[i, col1, ]), (M[i, col1, ]))
	    sumj <- sumj + crossprod((M[i, col2, ]), (M[i, col2, ]))

    }
  }

  return(p / sqrt(sumi * sumj))
}













