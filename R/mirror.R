#
# Code for mirroring
#   Two exported functions, mirrortrait and mirrorfill
#



#'Fill missing landmarks using mirrored values from other side of object
#'
#'Given an m x n x 3 matrix, replace a set of landmarks using their mirrored counterpart.
#'
#'@param a An m x n x 3 matrix containing 3D landmark data of n landmarks and m specimens
#'@param l1 Vector of indices for which landmarks to use to make a specimen midline
#'@param l2 Vector of indices for which landmarks to be 
#'
#'@details \code{l1} should be an even number length containing pairs of landmarks
#'  on either side of the specimen.
#'@export

mirrorfill <- function(a, l1, l2){
  stopifnot(is.numeric(a), dim(a)[3] == 3, length(dim(a)) == 3)

  # Count specimens and landmarks and check they're positive.
  m <- dim(a)[1]
  n <- dim(a)[2]
  stopifnot(m < 1, n < 1)

  # Make replicate that we will fill in
  a2 <- array(NA, dim = c(m, n, 3))

  # For each specimen, use mirrorfill1 to replace missing points.
  for(i in 1:m){
    a2[i, , ] <- mirrorfill1(a[i, , ], l1, l2)
  }

  return(a2)
}



#'Fill missing landmarks using mirrored values from other side of object
#'
#'Given an n x 3 matrix, replace a set of landmarks using their mirrored counterpark.
#'
#'@param s An n x 3 matrix containing 3D landmark data of n landmarks.
#'@param l1 Vector of indices for which landmarks to use to make a specimen midline.
#'@param l2 Vector of indices for which landmarks to be. 
#'
#'@details \code{l1} should be an even number length containing pairs of landmarks
#'  on either side of the specimen.
#'@export


mirrorfill1 <- function(s, l1, l2){
  # Check inputs
  stopifnot(is.numeric(a), dim(a)[2] == 3, length(dim(a)) == 2)

  # Count missing data points
  count1 <- sum(apply(a, 1, anyNA))

  # Initialise counter for number of landmarks replaced.
  count2 <- 0
  
  

}



mirrorfill1[s_,l1_,l2_] := Module[ {count1, count2, n, d, i, j, ii, jj},

    # counts missing data points
    count1 = validates[s];
    count2 = 0;
    {n,d} = midline[s, l1];
    If [!VectorQ[l2, IntegerQ],
	    Print["fatal error in mirrorfill[]: mirror list is bad"];
	    Abort[];
    ];

    ns = s;
    For [i = 1, i <= Length[l2], i++,
	    j = If[Mod[i,2]==1, i+1, i-1];
	    If [j > Length[l2],
	      Print["fatal error in mirrorfill[]: number of mirrored points is odd"];
	      Abort[];
	    ];

	    ii = l2[[i]];
	    jj = l2[[j]];
	    If [ii < 1 || ii > Length[s] || jj < 1 || jj > Length[s],
	      Print["fatal error in mirrorfill[]: mirror index is out of range"];
	      Abort[];
	    ];
	    If [ StringQ[ns[[ii]]], Continue[]];
	    If [!StringQ[ns[[jj]]], Continue[]];
	
	    ns[[jj]] = reflect[ ns[[ii]], n, d ];
	    count2++;
    ];

    Print["mirrorfill[] reconstructed ", count2, " out of ", count1, " missing landmarks"];
    Return[ns];
];



#'Find the midline of an object
#' 
#'@param s An n x 3 matrix of the 3D coordinates of n landmarks.
#'@param l1 Indices of the landmarks to use in creating the midline.

midline <- function(s, l1){
  
  # Perform some checks on inputs.
  stopifnot(is.numeric(s), dim(s)[2] == 3, length(dim(s)) == 2)
  stopifnot(is.integer(l1))

  if(any(l1 < 1) || any(l1 > length(dim(s)[1]))){
    stop("Midline index is out of range")
  }

  # take landmarks from l1 and remove missing data rows
  mlist <- s[l1, ]
  mlist <- mlist[complete.cases(mlist), ]

  # Check there's enough points to make a midline
  if(NROW(mlist) < 3) stop('Too much missing data to create midline')

  return(bestplane(mlist))
}


#'Calculate the least squares plane of some points
#'
#'@param l An n x 3 matrix of landmarks to be used in creating the plane
#'
#'@return A list containing a unit vector n and a scalar d.
#'  n.x = d is the least squares plane.

bestplane <- function(l){
  c <- lcentroid(l)
  nl <- lshift(l, -c)
  
  # Do the stuff that was in largestev function
  n <- eigen(-t(nl) %*% nl)$vectors[, 1]
  if(n[which.max(abs(n))] < 0){
    n <- -n
  }

  # Calculate how well the plane fits
  fit <- sum(sapply(1:NCOL(nl), function(i) n %*% nl[i, ]^2 ))
  message('Fit of midline plane is ', fit, ' with ', NCOL(nl), ' landmarks.')
  return(list(n, n %*% c))
}









