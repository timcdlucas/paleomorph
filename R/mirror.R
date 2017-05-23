#
# Code for mirroring
#   2 exported functions, mirrorfill1 and mirrorfill
#



#'Fill missing symmetrical landmarks for all specimens in an array using mirrored values from other side of a bilaterally symmetrical object where present
#'
#'Given an N x 3 x M matrix, where N is the number of landmarks, 3 is the number of dimensions, and M is the number of specimens, fill in missing landmarks using their mirrored counterpart.

#'
#'@param A An N x 3 x M matrix where N is the number of landmarks, 3 is the number of dimensions, and M is the number of specimens.
#'@param l1 Vector of indices for which landmarks to use to make a specimen midline
#'@param l2 Vector or matrix of pairs of symmetrical landmarks
#'
#'
#'@details \code{l2} should be either 
#'  \itemize{
#'    \item An even length vector containing pairs of landmarks on either side of the specimen. 
#'      i.e. l2[1] and l2[2] are paired, l2[3] and l2[4] are paired etc.
#'    \item A two column matrix with each row giving a pair of symmetrical landmarks.
#'  }
#'
#'@details \code{l2} should be an even number length containing pairs of landmarks
#'  on either side of the specimen.
#'@export
#'@examples
#' 
#' # Make objects that are in the z = 0 plane.
#' A <- array(rnorm(12 * 3 * 4), dim = c(12, 3, 4))
#' 
#' # Make it symmetrical
#' A[, 3, ] <- 0.1
#' A2 <- A
#' A2[, 3, ] <- -0.1
#' A <- abind::abind(A, A2, along = 1)
#' 
#' 
#' # Add some missing data points
#' missinga <- A
#' missinga <- abind::abind(missinga, array(NA, dim = c(2, 3, 4)), along = 1)
#' 
#' mirrorA <- mirrorfill(missinga, l1 = c(1:24), l2 = c(23, 25, 24, 26))
#' 


mirrorfill <- function(A, l1, l2){
  stopifnot(is.numeric(A), dim(A)[2] == 3, length(dim(A)) == 3)

  # Count specimens and landmarks and check they're positive.
  m <- dim(A)[3]
  n <- dim(A)[1]
  stopifnot(m > 1, n > 1)

  # Make replicate that we will fill in
  a2 <- array(NA, dim = c(n, 3, m))

  # For each specimen, use mirrorfill1 to replace missing points.
  for(i in 1:m){
    a2[, , i] <- mirrorfill1(A[, , i], l1, l2)
  }

  return(a2)
}



#'Fill missing landmarks for a single specimen using mirrored values from other side of object
#'
#'Given an n x 3 matrix, replace a set of landmarks using their mirrored counterpart.
#'
#'@param s An n x 3 matrix containing 3D landmark data of n landmarks.
#'@param l1 Vector of indices for which landmarks to use to make a specimen midline.
#'@param l2 Vector or matrix of pairs of symmetrical landmarks.
#'
#'
#'@details \code{l2} should be either 
#'  \itemize{
#'    \item An even length vector containing pairs of landmarks on either side of the specimen. 
#'      i.e. l2[1] and l2[2] are paired, l2[3] and l2[4] are paired etc.
#'    \item A two column matrix with each row giving a pair of symmetrical landmarks.
#'  }
#'@export
#'@examples
#' 
#' 
#' # Make data that is reflected in x plane
#' s <- matrix(c(
#'   rep(c(1, -1), each = 4),
#'   rep(c(4, 4, -4, -4), 2),
#'   rep(c(4, -4), 4)
#' ), ncol = 3)
#' 
#' 
#' # Now add some empty data
#' s <- rbind(s, NA)
#' 
#' # Mirror the NA point (row 9) using its complimentary landmark, point 5.
#' mirrorS <- mirrorfill1(s, l1 = 1:8, l2 = c(5, 9))




mirrorfill1 <- function(s, l1, l2){
  # Check inputs
  stopifnot(is.numeric(s), dim(s)[2] == 3, length(dim(s)) == 2)
  stopifnot(is.integer(l1), is.integer(l2) | is.numeric(l2))

  if(is.vector(l2)){
    if(length(l2) %% 2) stop('Number of mirrored points is odd')
  } else if(is.matrix(l2)){
    if(dim(l2)[2] != 2){
      stop('l2 should be a even length vector or a two column matrix.')
    }
  } else {
    stop('l2 should be a even length vector or a two column matrix.')
  } 

  # Get l2 into vector format
  #  I should probably rewrite whole function to use matrix form. Maybe later.
  if(is.matrix(l2)){
    l2 <- as.vector(t(l2))
  }

  # Count missing data points
  count1 <- sum(apply(s, 1, anyNA))

  # Initialise counter for number of landmarks replaced.
  count2 <- 0
  
  # Find the center plane
  mid <- midline(s, l1)
  
  # Make a copy of the shape
  ns <- s


  # for each element in the replacement list l2
  for(i in 1:length(l2)){
    # In list of integers a, b, c, d,
    #   a and b are a pair and c and d are a pair
    #   Compare each pair both ways and replace missing values with mirrored version
    #   i.e. if a has missing data, replace a with reflection of b
    #   but if b has missing data, replace b with reflection of a
    #   If neither has missing data do nothing. If both have missing data, do nothing.
    if(i %% 2 == 1){
      j <- i + 1
    } else {
     j <- i -1
    }
    # Get the indices for the landmarks    
    ii <- l2[i]
    jj <- l2[j]
    if(ii < 1 || ii > dim(s)[1] || jj < 1 || jj > dim(s)[1]){
      stop("fatal error in mirrorfill: mirror index is out of range")
    }

    # if landmark ii is complete data and jj has missing data, replace ii with mirror of jj.
    #   (the reverse replacement is tested in a later iteration of the loop)
    if(!anyNA(ns[ii, ]) & anyNA(ns[jj, ])){
      ns[jj, ] <- reflect(ns[ii, ], mid$n, mid$d)
      count2 <- count2 + 1
    }
  }
  print(paste("mirrorfill reconstructed ", count2, " out of ", count1, " missing landmarks"))
  return(ns)
  
}



#Find the midline of an object
# 
#@param s An n x 3 matrix of the 3D coordinates of n landmarks.
#@param l1 Indices of the landmarks to use in creating the midline.

midline <- function(s, l1){
  
  # Perform some checks on inputs.
  stopifnot(is.numeric(s), dim(s)[2] == 3, length(dim(s)) == 2)
  stopifnot(is.integer(l1))

  if(any(l1 < 1) || any(l1 > dim(s)[1])){
    stop("Midline index is out of range")
  }

  # take landmarks from l1 and remove missing data rows
  mlist <- s[l1, ]
  mlist <- mlist[complete.cases(mlist), ]

  # Check there's enough points to make a midline
  if(NROW(mlist) < 3) stop('Too much missing data to create midline')

  return(bestplane(mlist))
}
  

#Calculate the least squares plane of some points
#
#@param l An n x 3 matrix of landmarks to be used in creating the plane
#
#@return A list containing a unit vector n and a scalar d.
#  n.x = d is the least squares plane.

bestplane <- function(l){
  # centre the points
  c <- lcentroid(l)
  nl <- lshift(l, -c)
  
  n <- smallestev(eigen(cov(nl)))


  # Calculate how well the plane fits
  fit <- sum(sapply(1:NROW(nl), function(i) n %*% nl[i, ]^2 ))
  message('Fit of midline plane is ', fit, ' with ', NROW(nl), ' landmarks.')
  return(list(n = n, d = n %*% c))
}




# Reflects p in the plane defined by n.x = d
# n is a unit vector. d is a scalar.
# 
#@param p a length 3 vector (a 3D point in space)
#@param n Length three vector of coefficients for a plane nx = d
#@param d Length one vector giving coefficient d for a plane nx = d

reflect <- function(p, n, d){
  p - 2 * (n %*% p - d) * n/(n %*% n)
}




