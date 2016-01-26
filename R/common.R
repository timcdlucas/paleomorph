#' Decompose matric to a rotation and symetric positive definite matrix
#'
#' rpdecompose: Any invertible matrix can be uniquely decomposed
#'   M = RP, where R is a rotation and P is symmetric positive
#'   definite.  This returns the list {R,P}.
#'
#'
#'@param m A matrix
#'
#'@return A list with matrices R (the rotation) and P (the symetric 
#'  positive definite matrix) such that M = RP. 



rp_decompose <- function(m){
  stopifnot(is.numeric(m), !any(is.na(m)))
  if(nrow(m) != ncol(m)) stop('m must be a square matrix')

  sv <- svd(m)
  R <- t(sv$v) %*% sv$u
  rpdecompose <- list(R = R, P = t(R) %*% m)
  return(rpdecompose)
}


#' Common manipulations of 3D points
#'
#'  lshift: translates all 3D points in matrix A by a common displacement
#'
#'@param A An n x 3 matrix
#'@param v Length 3 displacement vector
#'
#'@return A matrix of the same dimensions as A
#'@name lshift

lshift <- function(A, v){
  A2 <- apply(A, 1, function(x)
          if(any(is.na(x))){
            x
          } else {
            x + v
          }
        )
  return(t(A2))
}



#'lrotate: rotates all 3D points in matrix A by a common matrix
#'
#'@param A An n x 3 matrix
#'@param m 3 x 3 rotation matrix
#'
#'@return A matrix of the same dimensions as A
#'@rdname lshift

lrotate <- function(A, m){
  A2 <- apply(A, 1, function(x)
          if(any(is.na(x))){
            x
          } else {
            x %*% m
          }
        )
  return(t(A2))  
}




#'lcentroid: computes centroid of all 3D points in matrix A (cannot have missing data)
#'
#'@param A An n x 3 matrix
#'
#'@return A matrix of the same dimensions as A
#'@rdname lshift

lcentroid <- function(A){
  stopifnot(!any(is.na(A)))
  A2 <- apply(A, 2, mean)
}



#'lcentroid2: like lcentroid but tolerates missing data
#'
#'@param A An n x 3 matrix
#'
#'@return A matrix of the same dimensions as A
#'@rdname lshift

lcentroid2 <- function(A){
  stopifnot(is.numeric(A), dim(A)[2] == 3, is.matrix(A))
  
  # Find mean for each spatial dimension after removing rows that have an NA.
  A2 <- apply(A[complete.cases(A), ], 2, mean)
  return(A2)
}




#'lnorm: returns total square length, tolerates missing data
#'
#'@param A An n x 3 matrix
#'
#'@return A scalar of total square length.
#'@rdname lshift

lnorm <- function(A){
  stopifnot(is.numeric(A), dim(A)[2] == 3, is.matrix(A))
  # rm missing data rows
  # Square each row
  # sum
    
  A2 <- sum(apply(A[complete.cases(A), ], 2, function(x) x %*% x))
  return(A2)
}




#'lscale: like lshift[], but multiplicative not additive
#'
#'@param A An n x 3 matrix
#'@param v Length 3 scaling vector
#'
#'@return A scalar of total square length.
#'@rdname lshift 

lscale <- function(A, v){
  A2 <- apply(A, 1, function(x)
          if(any(is.na(x))){
            x
          } else {
            x * v
          }
        )
  return(t(A2))
}

