% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSpecimens.R
\name{plotSpecimens}
\alias{plotSpecimens}
\title{Plot an array of specimen landmark data in an interactive 3D frame}
\usage{
plotSpecimens(A, l1 = NULL, midlineSpecimens = NULL, cols = NULL,
  bySpecimen = TRUE, planeOptions = NULL, ...)
}
\arguments{
\item{A}{An N x 3 x M array.}

\item{l1}{Optional vector of indices for which landmarks to use to make a specimen midline. If NULL, no midline plane is plotted.}

\item{midlineSpecimens}{Numeric vector indicating which specimens should be used to built the midline plane. If NULL, but l1 is defined, all specimens are used.}

\item{cols}{A vector of colours.}

\item{bySpecimen}{Logical that determined whether points should be coloured by specimen (default) or by landmark.}

\item{planeOptions}{Named list of parameters passed to \code{\link[rgl]{rgl.material}} to control the appearence of 
plotted mirror planes.}

\item{...}{Further parameters passed to \code{plot3d}.}
}
\description{
This function requires the rgl package.
  Given a N x 3 x M array (where M is the number of specimens
  and N is the number of landmarks), as used elsewhere in this package,
  plot each specimen in a different colour in an intereactive
  3D frame.
}
\examples{
A <- array(rep(rnorm(3 * 20, sd = 30), by = 6) + rnorm(6 * 20 * 3), 
       dim = c(20, 3, 6))
plotSpecimens(A)

plotSpecimens(A, bySpecimen = FALSE)

plotSpecimens(A, cols = grey(seq(0, 1, length.out = 6)))

plotSpecimens(A, l1 = c(1:4), planeOptions = list(alpha = 0.4, color = 'red'))


}
\seealso{
\code{\link[rgl]{plot3d}} \code{\link{mirrorfill}} \code{\link[rgl]{planes3d}} \code{\link[rgl]{rgl.material}}
}

