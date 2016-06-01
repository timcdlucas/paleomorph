
#' Count the number of missing landmarks in an array
#'
#'@param A An n x m x 3 array where n is the number of specimens and m is the number of landmarks.
#'
#'@return A length n vector giving the number of missing landmarks for each specimen.
#'@export
#'@examples
#'
#'
#'
#'  a <- array(1:(3*6*7), dim = c(7, 3, 6))
#'  a[2, , 1] <- NA
#'  countMissing(a)

countMissing <- function(A){
  stopifnot(is.numeric(A), length(dim(A)) == 3, dim(A)[2] == 3)
  completeLandmarks(A)
  
  # Find landmarks with missing data
  miss <- apply(A, c(3, 2), function(x) any(is.na(x)))

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


