context('deltaa')

test_that('deltaa works', {

  # delta should be 0
  olda <- array( rep(1:6, 3), dim = c(3, 2, 3))
  newa <- array( rep(1:6, 3), dim = c(3, 2, 3))

  m <- 3
  n <- 2

  expect_equal(0, deltaa(olda, newa, m, n))


  # small example.
  # each coord in newa is +1 coord in olda.
  # Each cooord has distance of 1.73205080756887719318.
  # Total is 1.732051 * 2 * 2

  olda <- array( rep(1:4, 3), dim = c(2, 2, 3))
  newa <- array( rep(2:5, 3), dim = c(2, 2, 3))

  m <- 2
  n <- 2

  expect_equal(1.73205080756887719318 * 2 * 2, deltaa(olda, newa, m, n))

})



test_that('deltaa works for negative coords', {

  # small example.
  # each coord in newa is +1 coord in olda.
  # Each cooord has distance of 1.73205080756887719318.
  # Total is 1.732051 * 2 * 2

  olda <- array( rep(0, 3), dim = c(2, 2, 3))
  newa <- array( rep(1, 3), dim = c(2, 2, 3))

  m <- 2
  n <- 2

  positive <- deltaa(olda, newa, m, n)

  newa <- array( rep(-1, 3), dim = c(2, 2, 3))
  negative <- deltaa(olda, newa, m, n) 
 
  expect_equal(positive, negative)


})


