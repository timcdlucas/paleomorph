
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
#' Given an M x N x 3 array with no missing data and a template
#'   M x N x 3 array with missing data, replace the missing data
#'   in the first array.
#'
#'@param a An M x N x 3 array.
#'@param b An M x N x 3 array to use a template. 
#'
#'@return An M x N x 3 array with NAs for missing data. 

unzapa <- function(a, b){

  a[is.na(b)] <- NA

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
#'@return An M x N x 3 array 


pcrstep <- function(a, maxiter = 1000, tolerance = 10e7){

  stopifnot(is.numeric(na))
  done <- FALSE

  # remove missing data
  na <- zapa(na)

  for(count in 1:maxiter){
    # Make copy
    na2 <- na
    if(count %% 50 == 0) message(paste('Iteration: ', count))

    for(i in 2:dim(na)[1]){
      # Compute ta, the temporary matrix that will be used
	    # in this iteration
      ta <- matrix(0, nrow = dim(na)[2], ncol = 3)
      for(j in 1:dim(na)[1]){
        if(i != j){
          ta = ta + na[j, , ] 
        }

 	    
	      # Compute na[i,,]^T (ta)
	      # The rotation which best approximates this matrix will be
	      # applied to na[i,,]
        c <- t(na[i, , ]) %*% ta

        # Take rotation part from Singular value decomposition
        r <- rp_decompose(c)$R
        # Apply rotation to na[i,,].
        na <- lrotate(na[i, , ], r)
      }
      # Test to see if approximation is good enough.
      #   Break out of loop if it is.
      if(deltaa(na2, na, dim(na2)[1], dim(na2)[2]) < tolerance) break()
        
    }     

  na <- unzapa(
  
    na = unzapa[na,a,m,n];
    Print["rstep: score=", scorea[na,m,n], " delta=", deltaa[a,na,m,n], " iterations=", count];
    Return[na];      
      
  }
  
}
  



(*
 * returns new "a"
 *)
pcrstep[a_,m_,n_] := Module[ {na, na2, count, ta, i, j, k, c, r},

    na = zapa[a,m,n];	(* no missing data in sight now *)

    done = False;
    For[count = 1, count <= 1000, count++,

	na2 = na;   (* save *)
	If[Mod[count, 50] == 0, Print["rstep: reached ", count, "th iteration"]];

	For[i = 2, i <= m, i++,

	    (*
	     * Compute ta, the temporary matrix that will be used
	     * in this iteration
	     *)
	    ta = Table[0.0, {j,1,n}, {k,1,3}];
	    For[j = 1, j <= m, j++,
		If[ i != j, ta += na[[j]] ];
	    ];

	    (*
	     * Compute na[[i]]^T (ta)
	     * The rotation which best approximates this matrix will be
	     * applied to na[[i]]
	     *)
	    c = Transpose[na[[i]]] . ta;
	    r = rpdecompose[c][[1]];
	    na[[i]] = lrotate[na[[i]], r];
	];

	If[deltaa[na,na2,m,n] < 10^(-7), Break[]];
    ];

    na = unzapa[na,a,m,n];
    Print["rstep: score=", scorea[na,m,n], " delta=", deltaa[a,na,m,n], " iterations=", count];
    Return[na];
];




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




