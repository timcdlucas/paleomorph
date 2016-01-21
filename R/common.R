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




