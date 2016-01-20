
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
      for(k in 1:k){
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




#rpdecompose <- function(m){
#  
#  
#}

#pcrstep[a_,m_,n_] 
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




