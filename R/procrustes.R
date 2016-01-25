
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
#'@param b An M x N x 3 array to use a template. 
#'
#'@return An M x N x 3 array with NAs for missing data. 

unzapa <- function(a, b){

  a[is.na(b)] <- NA

  return(a)
}



#' Rotate all shapes until optimally aligned.
#'
#' Given an M x N x 3 array (M = no of specimens, N = no of landmarks.) 
#'   this will find the optimal rotation for all shapes so that they are as
#'   aligned as possible.
#'
#'@param a An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@param maxiter Maximum number of iterations to attempt
#'@param tolerance Difference between two iterations that will cause the search to stop. 
#'
#'@return An M x N x 3 array of aligned shapes


pcrstep <- function(a, maxiter = 1000, tolerance = 10e-7){

  stopifnot(is.numeric(a))

  # remove missing data
  na <- zapa(a)

  for(count in 1:maxiter){
    # Make copy
    na2 <- na
    if(count %% 50 == 0) message(paste('Iteration: ', count))

    # For each specimen except fist
    #   rotate that specimen 
    for(i in 2:dim(na)[1]){
      # Compute ta, the temporary matrix that will be used
	    # in this iteration
      ta <- matrix(0, nrow = dim(na)[2], ncol = 3)
      
      # For each specimen
      for(j in 1:dim(na)[1]){
        if(i != j){
          ta <- ta + na[j, , ] 
        }
      }
    
      # Compute na[i,,]^T (ta)
      # The rotation which best approximates this matrix will be
      # applied to na[i,,]
      c <- t(na[i, , ]) %*% ta

      # Take rotation part from Singular value decomposition
      r <- rp_decompose(c)$R
      # Apply rotation to na[i,,].
      na[i, , ] <- lrotate(na[i, , ], r)
      
    }    

    # Does new rotations only change the matrix a tiny bit?
    if(deltaa(na2, na, dim(na2)[1], dim(na2)[2]) < tolerance) break()
  }
  # Replace missing values
  na <- unzapa(na, a)

  # print some output
  message(paste("rstep: score =", scorea(na, dim(na)[1], dim(na)[2]), ", delta =", deltaa(a, na, dim(na)[1], dim(na)[2]), ", iterations = ", count))
  return(na)
}
  




#' Translate all shapes until optimally aligned.
#'
#' Given an M x N x 3 array (M = no of specimens, N = no of landmarks.) 
#'   this will find the optimal translation for all shapes so that they are as
#'   aligned as possible.
#'
#'@param a An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'
#'@return An M x N x 3 array of aligned shapes


pctstep <- function(a){

  stopifnot(is.numeric(a), dim(a)[3] == 3, length(dim(a)) == 3)

#     * The linear algebra here is a little tricky.
#     * We have M equations in M unknowns, but one equation is redundant
#     * and the solution space is invariant under a common translation.
#     *
#     * For now, we deal with this by throwing away the first equation
#     * and forcing the first unknown to equal 0.
#     * A more robust solution might be: retain all M equations and solve by
#     * least squares, and force the sum of the unknowns to 0.

  # m numer of specimens, n number of landmarks
  m <- dim(a)[1]
  n <- dim(a)[2]

  # inialise objects
  c <- matrix(0, nrow = m - 1, ncol = m - 1)
  bx <- rep(0, m - 1) 
  by <- rep(0, m - 1) 
  bz <- rep(0, m - 1) 


  for(i in seq(2, m)){
    for(j in seq(m)){
      if (i != j) {
        for(k in seq(n)){
          if (all(!is.na(a[i, k, ])) && all(!is.na(a[j, k, ]))) {
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
  for(i in 2:dim(na)[1]){
    na[i, , ] <- lshift(a[i, , ], -v[i - 1, ])
  }
  message("tstep: score = ", scorea(na, dim(na)[1], dim(na)[2]), ", delta = ", deltaa(na, a, dim(na)[1], dim(na)[2]))
  return(na)
}






#' Resize all shapes until optimally aligned.
#'
#' Given an M x N x 3 array (M = no of specimens, N = no of landmarks.) 
#'   this will find the optimal resizing for all shapes so that they are as
#'   aligned as possible.
#'
#'@param a An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'
#'@return An M x N x 3 array of aligned shapes

#pcsstep <- function(a){

#  stopifnot(is.numeric(a), dim(a)[3] == 3, length(dim(a)) == 3)





#pcsstep[a_,m_,n_] := Module[ {na, d, v, i, j, k},

#    na = Table[
#	   lscale[ a[[i]],
#		   1/Sqrt[ lnorm[ lshift[ a[[i]], -lcentroid2[a[[i]]] ] ] ]
#	 ], {i,1,m}
#    ];

#    (*
#     * Compute D, the matrix whose largest eigenvalue will give the scalings
#     *)
#    d = Table[0.0, {i,1,m}, {j,1,m}];
#    For [i = 1, i <= m, i++,
#	For [j = 1, j <= m, j++,
#	    If [i == j, Continue[]];
#	    For [k = 1, k <= n, k++,
#		If [StringQ[na[[i,k]]] || StringQ[na[[j,k]]], Continue[]];
#		d[[i,i]] -= na[[i,k]] . na[[i,k]];
#		d[[i,j]] += na[[i,k]] . na[[j,k]];
#	    ];
#	];
#    ];
#    debugd = d;

#    (*
#     * The actual scaling is done here
#     *)
#    v = largestev[d];
#    debugv = v;
#    For[i = 1, i <= Length[v], i++,
#	If[v[[i]] > 0.0001, Continue[]];
#	Print["fatal error: procrustes scaling is negative!"];
#	Abort[];
#    ];
#    na = Table[lscale[na[[i]], v[[i]]], {i,1,m}];
#    Print["sstep: score=", scorea[na,m,n], " delta=", deltaa[a,na,m,n]];
#    Return[na];
#];



