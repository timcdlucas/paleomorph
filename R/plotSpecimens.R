

#' Plot an array of specimen landmark data in an interactive 3D frame
#'
#' This function requires the rgl package.
#'   Given a N x 3 x M array (where M is the number of specimens
#'   and N is the number of landmarks), as used elsewhere in this package,
#'   plot each specimen in a different colour in an intereactive
#'   3D frame.
#'
#'@param a An N x 3 x M array.
#'@param l1 Optional vector of indices for which landmarks to use to make a specimen midline. If NULL, no midline plane is plotted. 
#'@param midlineSpecimens Numeric vector indicating which specimens should be used to built the midline plane. If NULL, but l1 is defined, all specimens are used.
#'@param cols A vector of colours. 
#'@param bylandmark Logical that determined whether points should be coloured by specimen (default) or by landmark.
#'@param planeOptions Named list of parameters passed to \code{\link[rgl]{rgl.material}} to control the appearence of 
#'  plotted mirror planes.
#'@param ... Further parameters passed to \code{plot3d}. 
#'
#'@seealso \code{\link[rgl]{plot3d}} \code{\link{mirrorfill}} \code{\link[rgl]{planes3d}} \code{\link[rgl]{rgl.material}}
#'@export
#'
#'@examples
#' a <- array(rep(rnorm(3 * 20, sd = 30), by = 6) + rnorm(6 * 20 * 3), 
#'        dim = c(20, 3, 6))
#' plotSpecimens(a)
#'
#' plotSpecimens(a, bySpecimen = FALSE)
#'
#' plotSpecimens(a, cols = grey(seq(0, 1, length.out = 6)))
#'
#' plotSpecimens(a, l1 = c(1:4), planeOptions = list(alpha = 0.4, color = 'red'))
#'
#' 

plotSpecimens <- function(a, 
                          l1 = NULL, 
                          midlineSpecimens = NULL, 
                          cols = NULL, 
                          bySpecimen = TRUE, 
                          planeOptions = NULL,
                          ...) {

  # rgl can be a pain to install so it is in suggests, not imports.
  #   So if the user calls this function, need to check it is installed
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package rgl needed for this function to work. Please install it with install.packages('rgl').",
      call. = FALSE)
  }
  
  # Check array is correct form
  stopifnot(length(dim(a)) == 3, dim(a)[2] == 3, is.numeric(a), is.logical(bySpecimen))
  
  # Create midline if indices given
  if(!is.null(l1)){
    if(is.null(midlineSpecimens)) midlineSpecimens <- 1:dim(a)[3]

    splitArrayToList <- list()
    for(i in 1:length(midlineSpecimens)){
      splitArrayToList[[i]] <- a[l1, , midlineSpecimens[i]]
    }
    
    X <- do.call(abind::abind, list(splitArrayToList, along = 1))

    mid <- midline(X, 1:NROW(X))  
  }
  

  # Create colour palette if not given
  if(is.null(cols)){
    if(bySpecimen){
      cols <- 1:dim(a)[3]
    } else {
      cols <- 1:dim(a)[1]
    }
  }    

  # Sort out planeOptions
  planeParams <- makePlaneParamsList(mid, planeOptions)


  # Put the first specimen into vectors
  #  This is mostly a way to give reasonable axes labels, without blocking 
  #   xyzlabs from using ...
  x <- a[, 1, 1]
  y <- a[, 2, 1]
  z <- a[, 3, 1]

  # Plot with colours either bySpecimen or by landmark
  if(bySpecimen){
    # Do 3D plots
    rgl::plot3d(x, y, z, col = cols[1], ...) 

    for(i in seq_len(dim(a)[3])[-1]){
      rgl::plot3d(a[, 1, i], a[, 2, i], a[, 3, i], add = TRUE, col = cols[i], ...) 
    }
  } else {
    # Do 3D plots
    rgl::plot3d(x, y, z, col = cols, ...) 

    for(i in seq_len(dim(a)[3])[-1]){
      rgl::plot3d(a[, 1, i], a[, 2, i], a[, 3, i], add = TRUE, col = cols, ...) 
    }
  }

  # Plot a plane if asked for.
  if(!is.null(l1)){
    do.call(rgl::planes3d, planeParams)
  }
    
}




# Utility function for messing around with parameters for plane plotting.

makePlaneParamsList <- function(mid, planeOptions){

  planeParams <- list(a = mid$n, d = mid$d[1])
  if(is.null(planeOptions)){
    planeParams <- c(planeParams, alpha = 0.5)
  } else {
    # Check planeParams is named
    if(is.null(names(planeOptions))) stop('planeOptions should be a named list')

    # Hack. List of knowns params from ?rgl.material
    knownParameters <- c("color", "alpha", "lit", "ambient", "specular", "emission", "shininess", "smooth", "texture", "textype", "texmipmap", "texminfilter", "texmagfilter", "texenvmap", "front", "back", "size", "lwd", "fog", "point_antialias", "line_antialias", "depth_mask", "depth_test")

    # Any params in planeOptions that are not in rgl.material?
    unknownParams <- sapply(names(planeOptions), function(x) x %in% knownParameters)
    # Print warning, telling user which parameters are unknown
    if(!all(unknownParams)){
      warning("Unknown parameters in planeOptions. I don't recognise: ", paste(names(planeOptions)[!unknownParams], collapse = ', '))
    }
    # Combine planeOptions with mid ready for plotting.
    planeParams <- c(planeParams, planeOptions)
  }
}


