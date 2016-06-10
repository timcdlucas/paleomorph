# Decompose matric to a rotation and symetric positive definite matrix
#
# rpdecompose: Any invertible matrix can be uniquely decomposed
#   M = RP, where R is a rotation and P is symmetric positive
#   definite.  This returns the list {R,P}.
#
#
#@param m A matrix
#
#@return A list with matrices R (the rotation) and P (the symetric 
#  positive definite matrix) such that M = RP. 



rp_decompose <- function(m){
  stopifnot(is.numeric(m), !any(is.na(m)))
  if(nrow(m) != ncol(m)) stop('m must be a square matrix')

  sv <- svd(m)
  R <- sv$v %*% t(sv$u)
  rpdecompose <- list(R = R, P = t(R) %*% m)
  return(rpdecompose)
}


# Common manipulations of 3D points
#
#  lshift: translates all 3D points in matrix A by a common displacement
#
#@param A An n x 3 matrix
#@param v Length 3 displacement vector
#
#@return A matrix of the same dimensions as A
#@name lshift

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



#lrotate: rotates all 3D points in matrix A by a common matrix
#
#@param m 3 x 3 rotation matrix
#
#@return A matrix of the same dimensions as A
#@rdname lshift

# source code param documentation. As added rdname helpfiles, shouldn't duplicate.
#param A An n x 3 matrix
#param m 3 x 3 rotation matrix

lrotate <- function(A, m){
  A2 <- t(apply(A, 1, function(x) x %*% m))
  
  if(anyNA(A2)){
    A2[!stats::complete.cases(A)] <- A[!stats::complete.cases(A)]
  }
  
  return(A2)  
}


#lcentroid: computes centroid of all 3D points in matrix A (cannot have missing data)
#
#
#@return A matrix of the same dimensions as A
#@rdname lshift

# source code param documentation. As added rdname helpfiles, shouldn't duplicate.
#param A An n x 3 matrix

lcentroid <- function(A){
  stopifnot(!any(is.na(A)))
  A2 <- apply(A, 2, mean)
}



#lcentroid2: like lcentroid but tolerates missing data
#
#
#@return A matrix of the same dimensions as A
#@rdname lshift

# source code param documentation. As added rdname helpfiles, shouldn't duplicate.
#param A An n x 3 matrix

lcentroid2 <- function(A){
  stopifnot(is.numeric(A), dim(A)[2] == 3, is.matrix(A))
  
  # Find mean for each spatial dimension after removing rows that have an NA.
  A2 <- apply(A[stats::complete.cases(A), ], 2, mean)
  return(A2)
}




#lnorm: returns total square length, tolerates missing data
#
#
#@return A scalar of total square length.
#@rdname lshift

# source code param documentation. As added rdname helpfiles, shouldn't duplicate.
#param A An n x 3 matrix

lnorm <- function(A){
  stopifnot(is.numeric(A), dim(A)[2] == 3, is.matrix(A))
  # rm missing data rows
  # Square each row
  # sum
    
  A2 <- sum(apply(A[stats::complete.cases(A), ], 2, function(x) x %*% x))
  return(A2)
}




#lscale: like lshift[], but multiplicative not additive
#
#@return A scalar of total square length.
#@rdname lshift

# source code param documentation. As added rdname helpfiles, shouldn't duplicate.
#param A An n x 3 matrix
#param v A length 3 vector defining the scale factors.

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






#completeLandmarks: Check that all landmarks are either complete or all NA
#
#@param a An N x 3 x M array. M = no of specimens, N = no of landmarks.
#
#@return NULL


completeLandmarks <- function(a){
  # Check that all landmarks are either complete or all NA
  # Which landmarks x species have at least one na.
  mask <- apply(a, c(1, 3), function(x) anyNA(x))

  # Find all data for the landmark x species that have at least one na
  missingList <- which(mask, arr.ind = TRUE)
  allMissing <- matrix(NA, nrow = nrow(missingList), ncol = 3)
  
  if(nrow(missingList) > 0){
    for(r in 1:nrow(missingList)){
      allMissing[r, ] <- a[missingList[r, 1], , missingList[r, 2]]
    }
  }

  # If any of allMissing is not NA, then there is incomplete landmarks

  if(any(!is.na(allMissing))){
    partial <- which(apply(a, c(1, 3), function(x) anyNA(x) & !all(is.na(x))), arr.ind = TRUE)
    colnames(partial) <- c('specimen', 'landmark')
    message('Some landmarks are partially complete and partially missing. The landmarks are: ')
    print(partial)
    stop('Exiting due to partially missing landmarks.')
  }
  
}






# largestev
#(*
# * This might be a bit paranoid, but the Mathematica documentation doesn't
# * seem to guarantee that
# *    - the last eigenvalue returned is the largest
# *    - the eigenvectors are normalized
# *    - if all components of an eigenvector have the same sign, then the
# *      positive sign is the one chosen
# *)

largestev <- function(eig){
  n <- eig$vectors[, which.max(eig$values)]
  if(all(n < 0)) n <- -n
  return(n)
}


