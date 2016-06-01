

#' Plot an array of specimen landmark data in an interactive 3D frame
#'
#' This function requires the rgl package.
#'   Given a M x N x 3 array (where M is the number of speciments
#'   and N is the number of landmarks), as used elsewhere in this package,
#'   plot each specimen in a different colour in an intereactive
#'   3D frame.
#'
#'@param a An M x N x 3 array.
#'@param cols A vector of colours. 
#'@param bylandmark Logical that determined whether points should be coloured by specimen (default) or by landmark.
#'@param ... Further parameters passed to \code{plot3d}.
#'
#'@seealso \code{\link[rgl]{plot3d}}
#'@export
#'
#'@examples
#' a <- array(rep(rnorm(6 * 20, sd = 30), each = 6) + rnorm(6 * 20 * 3), 
#'        dim = c(6, 20, 3))
#' plotSpecimens(a)
#'
#' plotSpecimens(a, bySpecimen = FALSE)
#'
#' plotSpecimens(a, cols = grey(seq(0, 1, length.out = 6)))
#'
#' 

plotSpecimens <- function(a, cols = NULL, bySpecimen = TRUE, ...) {

  # rgl can be a pain to install so it is in suggests, not imports.
  #   So if the user calls this function, need to check it is installed
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package rgl needed for this function to work. Please install it with install.packages('rgl').",
      call. = FALSE)
  }
  
  # Check array is correct form
  stopifnot(length(dim(a)) == 3, dim(a)[3] == 3, is.numeric(a), is.logical(bySpecimen))
  
  if(is.null(cols)){
    if(bySpecimen){
      cols <- 1:dim(a)[1]
    } else {
      cols <- 1:dim(a)[2]
    }
  }    

  # Put the first specimen into vectors
  #  This is mostly a way to give reasonable axes labels, without blocking 
  #   xyzlabs from using ...
  x <- a[1, , 1]
  y <- a[1, , 2]
  z <- a[1, , 3]

  if(bySpecimen){
    # Do 3D plots
    rgl::plot3d(x, y, z, col = cols[1], ...) 

    for(i in 2:dim(a)[1]){
      rgl::plot3d(a[i, , 1], a[i, , 2], a[i, , 3], add = TRUE, col = cols[i], ...) 
    }
  } else {
        # Do 3D plots
    rgl::plot3d(x, y, z, col = cols, ...) 

    for(i in 2:dim(a)[1]){
      rgl::plot3d(a[i, , 1], a[i, , 2], a[i, , 3], add = TRUE, col = cols, ...) 
    }

  }
    
}
