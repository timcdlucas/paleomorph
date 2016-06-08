
#' Run procrustese analysis to align 3D shapes. 
#' 
#'@param a An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@param scale Logical indicating whether the size of the objects should be scales
#'  as well as rotated and translated.
#'@param maxiter Maximum number of iterations to attempt
#'@param tolerance Difference between two iterations that will cause the search to stop. 
#'@param scaleDelta Logical determining whether deltaa should be scaled by the total number of landmarks.
#'@export
#'@details
#' A number of computations are run until the difference between two iterations is less than \code{tolerance}.
#'   The more specimens and landmarks you have, the less each landmark is allowed to move before this tolerance
#'   is reached. Setting \code{scaleDelta = TRUE} will make the alignment run faster but have potentially less 
#'   well aligned results. But the alignment between a large and small array of shapes should be more comparable
#'   with \code{scaleDelta = TRUE}. However, preliminary tests imply that run time scales linearly with 
#'   \code{scaleDelta} set to \code{TRUE} or \code{FALSE}. 
#'
#'@return A new (N x 3 x M) array, where each 3d vector has been transformed
#'  in a per-specimen way.  The transformation is chosen to maximize,
#'  in the least-squares sense, the distances between specimens.
#'
#'@examples
#' # Make an array with 6 specimens and 20 landmarks
#' a <- array(rep(rnorm(6 * 20, sd = 20), each = 6) + rnorm(6 * 20 * 3), 
#'       dim = c(6, 20, 3))
#'
#' # Align the data (although it is already largely aligned)
#' aligned <- procrustes(a)
#' 
#' plotSpecimens(aligned)
#'
#' 
#'
#'

procrustes <- function(a, scale = TRUE, scaleDelta = FALSE, maxiter = 1000, tolerance = 10e-6){
  stopifnot(is.numeric(a), is.logical(scale), length(dim(a)) == 3, dim(a)[3] == 3)

  # Check that all landmarks are either complete or all NA
  completeLandmarks(a)
  
  na <- pcistep(a, scale)

  for(iter in 1:maxiter){
    # Save a copy of the array
    na2 <- na 

    # Do an iteration of procrustese
    na <- pctstep(na, scaleDelta)
    na <- pcrstep(na, tolerance = tolerance, scaleDelta = scaleDelta)
    if(scale) na <- pcsstep(na, scaleDelta)

    # Check progress
    if(deltaa(na2, na, dim(na2)[1], dim(na2)[2], scaleDelta) < tolerance) break()
  }

  # If it didn't converge give a warning
  if(iter == maxiter){
    warning('After ', maxiter, ' iterations the solution has not converged.')
    warning('Delta = ', round(deltaa(na2, na, dim(na2)[1], dim(na2)[2], scaleDelta), 4), ', tolerance = ', tolerance)
  }

  return(na)
}


#' Returns the sum of squares of distances that we're trying to minimize.
#' 
#'@param arr An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@param m No of specimens
#'@param n No of landmarks
#'
#'@return The sum of squares distance

scorea <- function(arr, m, n){
  stopifnot(is.numeric(arr), is.numeric(m), is.numeric(n))
  sumw <- 0
  # For each specimen (except last)
  for(i in 1:(m - 1)){
    # For the other specimens that haven't yet been compared
    for(j in (i + 1):m){

      # For each landmark
      for(k in 1:n){
          if(!anyNA(arr[k, , i]) & !anyNA(arr[k, , j])){
            w <- arr[k, , i] - arr[k, , j]
            sumw <- sumw + w %*% w
          }
      }
    }
  }
  return(drop(sumw))
}





#' Returns the sum of squares of the distances between "a1" and "a2".
#' 
#' For each landmark on each sample, find distance between location given
#'   in a1 and a2. 
#' Used to see when a1 and a2 are very similar. e.g. deltaa(olda, newa, 10, 20, scaleDelta = FALSE) < 10e-7
#'
#'@param olda An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@param newa An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@param m No of specimens
#'@param n No of landmarks
#'@param scaleDelta Logical determining whether deltaa should be scaled by the total number of landmarks.
#'@param zap Do the objects need to have missing data removed.
#'
#'@return The sum of squares distances (length 1 numeric) between all landmarks on all speciments.



deltaa <- function(olda, newa, m, n, scaleDelta, zap = TRUE){
  stopifnot(dim(newa) == dim(olda), is.numeric(newa), is.numeric(olda), is.numeric(m), 
            is.numeric(n), dim(newa) == c(n, 3, m), is.logical(scaleDelta))
  if(zap){
    olda <- zapa(olda)
    newa <- zapa(newa)
  }

  diff <- olda - newa
  
  if(scaleDelta){
    delta <- sum(apply(diff, c(1, 3), function(x) sqrt(x %*% x))) / (dim(olda)[1] * dim(olda)[3])
  } else {
    delta <- sum(apply(diff, c(1, 3), function(x) sqrt(x %*% x)))
  }
  return(delta)
}



#' Replaces missing data points with c(0, 0, 0)
#'
#' Given an N x 3 x M array, returns an N x 3 x M array with no missing data.
#'   M is the number of specimens and N is the number of landmarks.
#'
#'@param a An N x 3 x M array.
#'
#'@return An N x 3 x M array with no missing data. 

zapa <- function(a){
  stopifnot(is.numeric(a))

  w <- which(apply(a, c(1, 3), function(x) anyNA(x)), arr.ind = TRUE)
  a[w[, 1], , w[, 2]] <- 0

  return(a)
}


#' Puts missing data back in to a specimen x landmark array
#'
#' Given an N x 3 x M array, and a template defining which data were
#'   missing, returns an N x 3 x M array with NAs for missing data.
#'   M is the number of specimens and N is the number of landmarks.
#'
#'@param a An N x 3 x M array.
#'@param b An N x 3 x M array to use a template. 
#'
#'@return An N x 3 x M array with NAs for missing data. 

unzapa <- function(a, b){

  a[is.na(b)] <- NA

  return(a)
}



#' Rotate all shapes until optimally aligned.
#'
#' Given an N x 3 x M array (M = no of specimens, N = no of landmarks.) 
#'   this will find the optimal rotation for all shapes so that they are as
#'   aligned as possible.
#'
#'@param a An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@param maxiter Maximum number of iterations to attempt
#'@param tolerance Difference between two iterations that will cause the search to stop. 
#'@param scaleDelta Logical determining whether deltaa should be scaled by the total number of landmarks.
#'
#'@return An N x 3 x M array of aligned shapes


pcrstep <- function(a, maxiter = 1000, tolerance = 10e-7, scaleDelta){

  stopifnot(is.numeric(a))

  # remove missing data
  na <- zapa(a)

  for(count in 1:maxiter){
    # Make copy
    na2 <- na
    if(count %% 50 == 0) message('Iteration: ', count)

    # For each specimen except fist
    #   rotate that specimen 
    for(i in 2:dim(na)[3]){
      # Compute ta, the temporary matrix that will be used
	    # in this iteration
      ta <- matrix(0, nrow = dim(na)[1], ncol = 3)

      # Essentially calculate mean of other shapes
      for(j in 1:dim(na)[3]){
        if(i != j){ 
          ta <- ta + na[, , j] 
        }
      } 
    
      # Compute na[i,,]^T (ta)
      # The rotation which best approximates this matrix will be
      # applied to na[i,,]
      c <- base::crossprod(ta, na[, , i])
      # Same but slower. Highlights difference to comments in Anjalis code.
      # c <- t(ta) %*% na[i,,]

      # Take rotation part from Singular value decomposition
      r <- rp_decompose(c)$R
      # Apply rotation to na[i,,].
      na[, , i] <- lrotate(na[, , i], r)
      
    }    

    
    # print(deltaa(na2, na, dim(na2)[1], dim(na2)[2]))
    # Does new rotations only change the matrix a tiny bit?
    if(deltaa(na2, na, dim(na2)[3], dim(na2)[1], scaleDelta = scaleDelta, FALSE) < tolerance) break()
  }


  # print some output
  message("rstep: score = ", scorea(na, dim(na)[3], dim(na)[1]), 
    ", delta = ", deltaa(a, na, dim(na)[3], dim(na)[1], scaleDelta = scaleDelta), 
    ", iterations = ", count)

  # Replace missing values
  na <- unzapa(na, a)

  return(na)
}
  




#' Translate all shapes until optimally aligned.
#'
#' Given an N x 3 x M array (M = no of specimens, N = no of landmarks.) 
#'   this will find the optimal translation for all shapes so that they are as
#'   aligned as possible.
#'
#'@param a An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@param scaleDelta Logical determining whether deltaa should be scaled by the total number of landmarks.
#'
#'@return An N x 3 x M array of aligned shapes


pctstep <- function(a, scaleDelta){

  stopifnot(is.numeric(a), dim(a)[2] == 3, length(dim(a)) == 3)

#     * The linear algebra here is a little tricky.
#     * We have M equations in M unknowns, but one equation is redundant
#     * and the solution space is invariant under a common translation.
#     *
#     * For now, we deal with this by throwing away the first equation
#     * and forcing the first unknown to equal 0.
#     * A more robust solution might be: retain all M equations and solve by
#     * least squares, and force the sum of the unknowns to 0.

  # m numer of specimens, n number of landmarks
  m <- dim(a)[3]
  n <- dim(a)[1]

  # inialise objects
  c <- matrix(0, nrow = m - 1, ncol = m - 1)
  bx <- rep(0, m - 1) 
  by <- rep(0, m - 1) 
  bz <- rep(0, m - 1) 


  for(i in seq(2, m)){
    for(j in seq(m)){
      if (i != j) {
        for(k in seq(n)){
          if (all(!is.na(a[k, , i])) && all(!is.na(a[k, , j]))) {
            c[i - 1, i - 1] <- c[i - 1, i - 1] + 1
            if (j > 1) {
              c[i - 1, j - 1] <- c[i - 1, j - 1] - 1
            }

            bx[i - 1] <- bx[i - 1] + a[i, k, 1] - a[j, k, 1]
            by[i - 1] <- by[i - 1] + a[i, k, 2] - a[j, k, 2]
            bz[i - 1] <- bz[i - 1] + a[i, k, 3] - a[j, k, 3]
          }          
        }
      }
    }
  }

	v <- t(rbind(solve(c, bx), solve(c, by), solve(c, bz)))

  na <- a
  for(i in 2:dim(na)[3]){
    na[, , i] <- lshift(a[, , i], -v[i - 1, ])
  }
  message("tstep: score = ", scorea(na, dim(na)[3], dim(na)[1]), ", delta = ", deltaa(na, a, dim(na)[3], dim(na)[1], scaleDelta))
  return(na)
}






#' Resize all shapes until optimally aligned.
#'
#' Given an N x 3 x M array (M = no of specimens, N = no of landmarks.) 
#'   this will find the optimal resizing for all shapes so that they are as
#'   aligned as possible.
#'
#'@param a An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@param scaleDelta Logical determining whether deltaa should be scaled by the total number of landmarks.
#'
#'
#'@return An N x 3 x M array of aligned shapes

pcsstep <- function(a, scaleDelta){

  stopifnot(is.numeric(a), dim(a)[2] == 3, length(dim(a)) == 3)

  # m numer of specimens, n number of landmarks
  m <- dim(a)[3]
  n <- dim(a)[1]

  na <- array(NA, dim = dim(a))

  for(i in 1:m){
    na[, , i] <- lscale(a[, , i], 1/sqrt(lnorm(lshift(a[, , i], -lcentroid2(a[, , i])))))
  }

  # Compute D, the matrix whose largest eigenvalue will give the scalings
  d <- matrix(0, nrow = m, ncol = m)
  # For each element of d (i.e. m x m)
  for(i in 1:m){
    for(j in 1:m){
      # ignore the diagonal
      if(j != i){
        # Loop through landmarks 
        for(k in 1:n){
          # ignore if the landmark in either specimens has missing data.
          if(!anyNA(na[k, , i]) && !anyNA(na[k, , j])){
            d[i, i] <- d[i, i] - (na[k, , i] %*% na[k, , i])
            d[i, j] <- d[i, j] + (na[k, , i] %*% na[k, , j])
          }
        }
      }  
    }
  }


  # Do the stuff that was in largestev function
  v <- eigen(d)$vectors[, 1]
  if(v[which.max(abs(v))] < 0){
    v <- -v
  }

  # The actual scaling is done here
  # Check that scaling isn't negative  
  if(any(v <= 0)){
    stop('Procrustes scaling is negative!')
  }

  for(i in 1:m){
    na[, , i] <- lscale(na[, , i], v[i])
  }
  
  message("sstep: score = ", scorea(na, m, n), ", delta = ", deltaa(a, na, m, n, scaleDelta))
  return(na)
}





#'Shifts each centroid to the origin.  This is not guaranteed
#'  to decrease the value of the objective function, so it makes
#'  no sense in later iterations; we just do it initially to get
#'  a head start on the convergence.
#'
#'  We also normalize the centroid size of each specimen to 1/Sqrt[m],
#'  so that the "total size" constraint is satisfied and scorea[] will
#'  become meaningful.
#'
#' Given an N x 3 x M array (M = no of specimens, N = no of landmarks.) 
#'   this will find the optimal resizing for all shapes so that they are as
#'   aligned as possible.
#'
#'@param a An N x 3 x M array. M = no of specimens, N = no of landmarks.
#'@param scale Logical indicating whether the size of the objects should be scaled.
#'
#'@return An N x 3 x M array of aligned shapes

pcistep <- function(a, scale = TRUE){

  na <- a

  # Shift centroid to origin
  for(i in 1:dim(na)[1]){
    na[i, , ] <- lshift(na[i, , ], -lcentroid2(na[i, , ]))
  }
  
  # Scale size of each specimen to 1/sqrt(m)
  if(scale){
    for(i in 1:dim(na)[1]){
      na[i, , ] <- lscale(na[i, , ],  1/sqrt(dim(na)[1] * lnorm(na[i, , ])))
    }
  }


  message("istep: score = ", scorea(zapa(na), dim(na)[1], dim(na)[2]))

  return(na)
}









