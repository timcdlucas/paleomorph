context('Test all functions used for mirroring')



test_that('Best plane works', {

  #a <- array(1:(3*40), dim = c(4, 10, 3))
  s <- matrix(1:21, ncol = 3, nrow = 7)

  m <- midline(s, 1:4)
  
  expect_equal(sum(m$n^2), 1)
  expect_true(length(m$n) == 3)
  expect_true(length(m$d) == 1)
  expect_true(all(names(m) == c('n', 'd')))
  expect_true(all( sapply(m, is.numeric)))

})


test_that('Best plane returns the best plane.', {

  # Make an object that is reflected in the horizontal, z = 0 plane.
  s <- matrix(rep(1:21, 2), byrow = TRUE, ncol = 3)
  s[1:7, 3] <- -s[1:7, 3]

  #s[, 3] <- s[, 3] + 5
  

  pl <- bestplane(s)
  
  expect_equal(drop(pl$n[1] + pl$n[2] - pl$d), 0.0)

})

test_that('midline returns best plane', {


})



test_that('bestplane and midline return same values when no missing data.', {


  # Make an object that is reflected in the horizontal, z = 0 plane.
  s <- matrix(rep(1:21, 2), byrow = TRUE, ncol = 3)
  s[1:7, 3] <- -s[1:7, 3]

  p1 <- bestplane(s)
  l1 <- 1:14
  p2 <- midline(s, l1)

  expect_equal(p1, p2)




  # Make an object that is reflected in the horizontal, z = 0 plane.
  s <- matrix(rep(rnorm(21), 2), byrow = TRUE, ncol = 3)
  s[1:7, 3] <- -s[1:7, 3]

  p1 <- bestplane(s)
  l1 <- 1:14
  p2 <- midline(s, l1)

  expect_equal(p1, p2)


})

test_that('Mirrorfill1 replaces points correctly.', {
  
  # Make an object that is reflected in the horizontal, z = 0 plane.
  s <- matrix(rep(1:21, 2), byrow = TRUE, ncol = 3)
  s[1:7, 3] <- -s[1:7, 3]

  # Now make a copy with missing data
  sna <- s
  sna[1, ] <- NA

  mirrorS <- mirrorfill1(sna, l1 = c(2:7, 9:14), l2 = c(1, 8))


})


test_that('Reflect works', {

  # Define a horizontal plane through origin
  n <- c(0, 0, 1)
  d <- 0
  
  # Test point (1, 1, 1). Should reflect to (1, 1, -1).
  p <- c(1, 1, 1)
  
  p2 <- reflect(p, n, d)
  
  expect_equal(p2, c(1, 1, -1))

})

