
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

dist <- function(v1,v2){
  return(sqrt((v1-v2)%*%(v1-v2)))
}


deltaa = function(a,na,m,n){
  delta = 0
  for(i in seq(m)){
    for(j in seq(n)){
      if (sum(a[i,j,]) != 0 && sum(na[i,j,]) != 0){
        delta = delta + dist(a[i,j,],na[i,j,])
      }
    }
  }
  return(delta)
}


rpdecompose <- function(m){
  
  
}

pcrstep[a_,m_,n_] 
pcrstep = function(na,m,n){
  for(count in seq(1000)){
    na2 = na
    for(i in seq(2,m)){
      ta = matrix(0,nrow=n,ncol=3)
      for(j in seq(m)){
        if(i != j){
          ta = ta + na[j,,] 
        }
        c = t(na[i,,]) %*% ta
        r = rpdecompose
      }
    }
    
  }
  
}



pctstep <- function(a, m, n){
  for(i in seq(2, m)){
    for(j in seq(m)){
      if (i != j) {
        for(k in seq(n)){
          if (sum(a[i, k, ]) != 0 && sum(a[j, k, ]) != 0) {
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
}




