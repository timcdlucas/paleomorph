
#' Count the number of missing landmarks in an array
#'
#'@param A An n x m x 3 array where n is the number of specimens and m is the number of landmarks.
#'
#'@return A length n vector giving the number of missing landmarks for each specimen.
#'@export

countMissing <- function(A){
  stopifnot(is.numeric(A), length(dim(A)) == 3, dim(A)[3] == 3)
  
  # Find landmarks with missing data
  miss <- apply(A, c(1, 2), function(x) any(is.na(x)))

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
  
  counts
  
}


