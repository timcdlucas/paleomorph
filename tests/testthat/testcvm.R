context('Test that all cvm functions work.')

test_that('dotcorr calculates congruence coefficient properly.', {

M <- array(0, dim = c(4, 2, 3))

# Landmarks are perfectly correlated
M[, 1, ] <- rep(1:4, by = 3)
M[, 2, ] <- rep(1:4, by = 3)




M <- array(0, dim = c(4, 2, 3))

# Landmarks are perfectly correlated
M[, 1, ] <- 1:12
M[, 2, ] <- 1:12

r <- dotcorrentry(M, 1, 2)

#expect_equal(r, 1)

})


test_that('Congruence is never outside -1, 1', {

  set.seed(999)
  r <- replicate(100, {
    # I don't think they need to be centered first. So use points not around (0, 0, 0)
    mu <- rnorm(1, 0, 3)
    M <- array(rnorm(4*2*3, mu), dim = c(4, 2, 3))
    dotcorr(M)
  })

  expect_true(all(r >= -1))
  expect_true(all(r <= 1))



})

