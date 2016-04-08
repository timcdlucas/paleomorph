context('Test the utility functions.')



test_that('countMissing counts and prints correct information.', {
  
  set.seed(2)

  A <- array(rnorm(10 * 20 * 3), dim = c(10, 20, 3))
  
  
  A[1, c(3, 4, 5), ] <- NA
  A[5, 8, ] <- NA

  miss <- countMissing(A)

  expect_equal(miss, c(3, 0, 0, 0, 1, rep(0, 5)))

})

