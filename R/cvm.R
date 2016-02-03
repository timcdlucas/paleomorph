
#' Calculate 1D covariance matrix 
#' 
#' Each row of M should correspond to one _specimen _.  Each column of
#' M should correspond to one _landmark _.  In the 1D case, each element
#' of M should be either a number or a string such as "?" which indicates
#' missing data.  In the 3D case, each element of M should be either
#' a 3-component vector or a string such as "?" which indicates missing
#' data.
#'
#'@param M An M x N matrix. M = no of specimens, N = no of landmarks.
#'@export
#'
#'@return 1D covariance matrix





#' Calculate 3D covariance matrix 
#' 
#' Each row of M should correspond to one _specimen _.  Each column of
#' M should correspond to one _landmark _.  In the 1D case, each element
#' of M should be either a number or a string such as "?" which indicates
#' missing data.  In the 3D case, each element of M should be either
#' a 3-component vector or a string such as "?" which indicates missing
#' data.
#'
#'@param M An M x N x 3 array. M = no of specimens, N = no of landmarks.
#'@export
#'
#'@return 3D covariance matrix

dotcvm <- function(M){
  N <- NA


}



dotcvm[M_] := Module[ {i,j,N,e},
   N = Table[dotcvmentry[M,i,j], {i,1,Dimensions[M][[2]]}, {j,1,Dimensions[M][[2]]}];
   e = Min[Eigenvalues[N]];
   If[e < 0, Print["warning: CVM has negative eigenvalue ", e];];
   N
];


# Check columns have enough data and then calculate covariance between columns
dotcvmentry <- function(M, col1, col2){
  n <- 0
  s1 <- c(0, 0, 0)
  s2 <- c(0, 0, 0)

  for(i in 1:dim(M)[2]){
    if(!anyNA(M[, c(col1, col2), ])){
      print(i)
      n <- n + 1
      s1 <- s1 + M[i, col1, ]
      s2 <- s2 + M[i, col1, ]
    }
  }

  if(n < 1) stop(paste("There is too much missing data  covary columns", col1, "and", col2))

  s1 <- s1/n
  s2 <- s2/n

  p <- 0
  for(i in 1:dim(M)[1]){
    if(!anyNA(M[, c(col1, col2)])){
      p <- p + crossprod((M[i, col1, ] - s1), (M[i, colw, ] - sw))
    }
  }

  return(p/(n - 1))
}

dotcvmentry[M_, col1_, col2_] := Module[ {i,n,s1,s2,p},
   n = 0;
   s1 = {0,0,0};
   s2 = {0,0,0};
   For[i = 1, i <= Dimensions[M][[1]], i++,
      If[ !StringQ[M[[i,col1]]] && !StringQ[M[[i,col2]]],
         n++;
         s1 = s1 + M[[i,col1]];
         s2 = s2 + M[[i,col2]];
      ]
   ];
   If[n<=1,
      Print["there is too much missing data to covary columns ", col1, " and ", col2];
      Abort[];
   ];
   s1 = s1/n;
   s2 = s2/n;

   p = 0;
   For[i = 1, i <= Dimensions[M][[1]], i++,
      If[ !StringQ[M[[i,col1]]] && !StringQ[M[[i,col2]]],
         p = p + (M[[i,col1]] - s1) . (M[[i,col2]] - s2);
      ];
   ];
   p/(n-1)
];






