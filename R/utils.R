
#' Count the number of missing landmarks in an array
#'
#'@param A An N x 3 x M array where N is the number of landmarks, 3 is the number of dimensions, and M is the number of specimens.
#'
#'@return A length n vector giving the number of missing landmarks for each specimen.
#'@export
#'@examples
#'
#'
#'
#'  A <- array(1:(3*6*7), dim = c(7, 3, 6))
#'  A[2, , 1] <- NA
#'  countMissing(A)

countMissing <- function(A){
  stopifnot(is.numeric(A), length(dim(A)) == 3, dim(A)[2] == 3)
  completeLandmarks(A)
  
  # Find landmarks with missing data
  miss <- apply(A, c(3, 1), function(x) any(is.na(x)))

  # Get counts for each specimen
  counts <- rowSums(miss)

  # How many at least one missing
  nIncomplete <- sum(counts > 0)

  # MAX
  whichMax <- which.max(counts)
  maxMissing <- max(counts)

  message('')
  message(paste(nIncomplete, 'specimens have missing data.'))
  message(paste('Max missing: specimen', whichMax, 'with', maxMissing, 'missing landmarks.\n'))
  
  return(counts)
  
}


