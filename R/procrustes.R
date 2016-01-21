
#' Returns the sum of squares of distances that we're trying to minimize.
#' 
#'@param mat An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@param m No of specimens
#'@param n No of landmarks
#'
#'@return The sum of squares distance

scorea <- function(arr, m, n){
  stopifnot(is.numeric(arr), !is.na(arr), is.numeric(m), is.numeric(n))
  sumw <- 0
  # For each specimen (except last)
  for(i in 1:(m - 1)){
    # For the other specimens that haven't yet been compared
    for(j in (i + 1):m){

      # For each landmark
      for(k in 1:n){
          w <- arr[i, k, ] - arr[j, k, ]
          sumw <- sumw + w %*% w
      }
    }
  }
  return(drop(sumw))
}





#' Returns the sum of squares of the distances between "a1" and "a2".
#' 
#' For each landmark on each sample, find distance between location given
#'   in a1 and a2. 
#' Used to see when a1 and a2 are very similar. e.g. deltaa(olda, newa, 10, 20) < 10e-7
#'
#'@param olda An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@param newa An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@param m No of specimens
#'@param n No of landmarks
#'
#'@return The sum of squares distances (length 1 numeric) between all landmarks on all speciments.


deltaa <- function(olda, newa, m, n){
  stopifnot(dim(newa) == dim(olda), is.numeric(newa), is.numeric(olda), is.numeric(m), is.numeric(n), dim(newa) == c(m, n, 3))
  delta <- 0
  for(i in seq(m)){
    for(j in seq(n)){
      delta <- delta + stats::dist(rbind(olda[i, j, ], newa[i, j, ]))
    }
  }
  return(as.vector(delta))
}


#' Replaces missing data points with c(0, 0, 0)
#'
#' Given an M x N x 3 array, returns an M x N x 3 array with no missing data.
#'   M is the number of specimens and N is the number of landmarks.
#'
#'@param a An M x N x 3 array.
#'
#'@return An M x N x 3 array with no missing data. 

zapa <- function(a){
  stopifnot(is.numeric(a))
  for(i in 1:dim(a)[1]){
    for(j in 1:dim(a)[2]){
      if(any(is.na(a[i, j, ]))){
        a[i, j, ] <- c(0, 0, 0)
      }
    }
  }
  return(a)
}



#' Puts missing data back in to a specimen x landmark array
#'
#' Given an M x N x 3 array, and a template defining which data were
#'   missing, returns an M x N x 3 array with NAs for missing data.
#'   M is the number of specimens and N is the number of landmarks.
#'
#'@param a An M x N x 3 array.
#'@param b An M x N logical matrix with TRUEs where the data were missing.
#'
#'@return An M x N x 3 array with NAs for missing data. 

unzapa <- function(a, b){

  # Check that the template b and the array a match.
  matches <- TRUE
  for(i in 1:nrow(b)){
    for(j in 1:ncol(b)){
      if(any(a[i, j, ] != 0) & b[i, j]){
        # If there are any non zeroes where we think the data is missing, 
        #   set flag to false and break out of loop.
        matches <- FALSE
        break()
      }
      # If there are mismatches, break out of outer loop as well.
      if(!matches) break()
    }
  }  
  # Give a warning if there are mismatches.
  if(!matches){
    warning("Non-zeros in positions marked as missing. \nAre you sure they're missing?")
  }

  # loop over specimens.
  for(i in 1:nrow(b)){
    # loop over landmarks
    for(j in 1:ncol(b)){
      # For every specimen, landmark position in b, if it's true (i.e. marking missing)
      #   set all three dimension in a as NA.
      if(b[i, j]){
        a[i, j, ] <- c(NA, NA, NA)
      }
    }
  }
  return(a)
}


#pcrstep = function(na,m,n){
#  for(count in seq(1000)){
#    na2 = na
#    for(i in seq(2,m)){
#      ta = matrix(0,nrow=n,ncol=3)
#      for(j in seq(m)){
#        if(i != j){
#          ta = ta + na[j,,] 
#        }
#        c = t(na[i,,]) %*% ta
#        r = rpdecompose
#      }
#    }
#    
#  }
#  
#}



#pctstep <- function(a, m, n){
#  for(i in seq(2, m)){
#    for(j in seq(m)){
#      if (i != j) {
#        for(k in seq(n)){
#          if (sum(a[i, k, ]) != 0 && sum(a[j, k, ]) != 0) {
#            c[i - 1, i - 1] <- c[i - 1, i - 1] + 1
#            if (j > 1) {
#              c[i - 1, j - 1] <- c[i - 1, j - 1] - 1
#            }

#            bx[i - 1] <- bx[i - 1] + a[i, k, 1] - a[j, k, 1]
#            by[i - 1] <- by[i - 1] + a[i, k, 2] - a[j, k, 2]
#            bz[i - 1] <- bz[i - 1] + a[i, k, 3] - a[j, k, 3]
#          }          
#        }
#      }
#    }
#  }
#}




