context('Test all functions used for mirroring')


set.seed(991)
test_that('Best plane works', {


  s <- matrix(1:21, ncol = 3, nrow = 7)

  m <- midline(s, 1:4)
  
  expect_equal(sum(m$n^2), 1)
  expect_true(length(m$n) == 3)
  expect_true(length(m$d) == 1)
  expect_true(all(names(m) == c('n', 'd')))
  expect_true(all( sapply(m, is.numeric)))


})


test_that('Best plane returns the best plane.', {

  set.seed(3)
  # Make an object that is in the horizontal, z = 0 plane.
  s <- cbind(rnorm(10), rnorm(10), 0)

  pl <- bestplane(s)

  expect_equal(pl$n, c(0,0,1))  
  expect_equal(drop(pl$n[1] + pl$n[2] - pl$d), 0.0)


})

test_that('midline returns best plane', {

  # Make an object that is in the horizontal, z = 0 plane.
  s <- cbind(rnorm(14), rnorm(14), 0)

  pl <- midline(s, 1:14)
  pl2 <- midline(s, c(2:5, 9:12))

  expect_equal(drop(pl$n[1] + pl$n[2] - pl$d), 0.0)
  expect_equal(drop(pl2$n[1] + pl2$n[2] - pl2$d), 0.0)

})



test_that('bestplane and midline return same values when no missing data.', {


  # Make an object that is in the horizontal, z = 0 plane.
  s <- cbind(rnorm(14), rnorm(14), 0)


  p1 <- bestplane(s)
  l1 <- 1:14
  p2 <- midline(s, l1)

  expect_equal(p1, p2)



  set.seed(99)
  # Make an object that is in the horizontal, z = 0 plane.
  s <- cbind(rnorm(14), rnorm(14), 0)

  p1 <- bestplane(s)
  l1 <- 1:14
  p2 <- midline(s, l1)

  expect_equal(p1, p2)



})

test_that('Mirrorfill1 replaces points correctly.', {
  
  # Make an object that is in the horizontal, z = 0 plane.
  s <- cbind(rnorm(14), rnorm(14), 0)

  # Now add a data point not on the z = 0 plane
  s <- rbind(s, c(1, 2, 1), c(NA, NA, NA))

  mirrorS <- mirrorfill1(s, l1 = 1:14, l2 = c(15, 16))

  expect_equal(mirrorS[16, ], c(1, 2, -1))

  
  
  
  # Make an object that is in the horizontal, z = 0 plane.
  s <- cbind(runif(14), runif(14), c(-0.1, 0.1))
  sneg <- s
  sneg[, 1:2] <- -sneg[, 1:2]
  s <- rbind(s, sneg)
  
  # Now add a data point not on the z = 0 plane
  s <- rbind(s, c(1, 2, 1), c(NA, NA, NA))
  
  mirrorS <- mirrorfill1(s, l1 = 1:28, l2 = c(29, 30))
  
  expect_equal(mirrorS[30, ], c(1, 2, -1))
  
})



test_that('Mirrorfill1 replaces points correctly in another plane.', {
  
  
  # Make an object that is in the horizontal, y = 0 plane.
  s <- cbind(rnorm(14), 0, rnorm(14))
  
  # Now add a data point not on the z = 0 plane
  s <- rbind(s, c(1, 2, 1), c(NA, NA, NA))
  
  mirrorS <- mirrorfill1(s, l1 = 1:14, l2 = c(15, 16))
  
  expect_equal(mirrorS[16, ], c(1, -2, 1))
})


test_that('Mirrorfill1 replaces points correctly in third plane.', {
  
  
  # Make an object that is in the horizontal, y = 0 plane.
  s <- cbind(0, rnorm(14), rnorm(14))
  
  # Now add a data point not on the z = 0 plane
  s <- rbind(s, c(1, 2, 1), c(NA, NA, NA))
  
  mirrorS <- mirrorfill1(s, l1 = 1:14, l2 = c(15, 16))
  
  expect_equal(mirrorS[16, ], c(-1, 2, 1))
})



test_that('Mirrorfill1 replaces points approximately correctly', {
  
  
  # Make an object that is in the horizontal, y = 0 plane.
  s <- cbind(0.1, rnorm(30), rnorm(30))
  s2 <- s
  s2[, 1] <- - s2[, 1]
  s <- rbind(s, s2)
  
  # Now add a data point not on the x = 0 plane
  s <- rbind(s, c(1, 2, 1), c(NA, NA, NA))
  
  mirrorS <- mirrorfill1(s, l1 = 1:60, l2 = c(61, 62))
  
  expect_equal(mirrorS[62, ], c(-1, 2, 1))
})


test_that('Example from mirrorfill1 works correctly', {
  
  # Make data that is reflected in x plane
  s <- matrix(c(
    rep(c(1, -1), each = 4),
    rep(c(4, 4, -4, -4), 2),
    rep(c(4, -4), 4)
  ), ncol = 3)

  
  # Now add some empty data
  s <- rbind(s, NA)
  
  # We will pair the NA row with row 5, which mirrors row 1. Row 1 is our test row.
  stest <- s[1, ]
  
  
  # Mirror point 1 using it's complimentary landmark, point 8.
  mirrorS <- mirrorfill1(s, l1 = 1:8, l2 = c(5, 9))
  expect_equal(mirrorS[9, ], stest)
  
})


test_that('Mirrorfill1 replaces points correctly with l1 as 2 col matrix.', {
  
  # Make an object that is in the horizontal, z = 0 plane.
  s <- cbind(rnorm(14), rnorm(14), 0)

  # Now add a data point not on the z = 0 plane
  s <- rbind(s, c(1, 2, 1), c(NA, NA, NA), c(3, 4, -2), c(NA, NA, NA))

  mirrorS <- mirrorfill1(s, l1 = 1:14, l2 = matrix(15:18, byrow = TRUE, ncol = 2))

  expect_equal(mirrorS[16, ], c(1, 2, -1))
  expect_equal(mirrorS[18, ], c(3, 4, 2))

})


test_that('Reflect works', {

  # Define a horizontal plane through origin
  n <- c(0, 0, 1)
  d <- 0
  
  # Test point (1, 1, 1). Should reflect to (1, 1, -1).
  p <- c(1, 1, 1)
  
  p2 <- reflect(p, n, d)
  
  expect_equal(p2, c(1, 1, -1))

  
  n <- c(0, 1, 0)
  d <- 0
  
  # Test point (1, 1, 1). Should reflect to (1, 1, -1).
  p <- c(1, 1, 1)
  
  p2 <- reflect(p, n, d)
  
  expect_equal(p2, c(1, -1, 1))
  
  
  n <- c(1, 0, 0)
  d <- 0
  
  # Test point (1, 1, 1). Should reflect to (1, 1, -1).
  p <- c(1, 1, 1)
  
  p2 <- reflect(p, n, d)
  
  expect_equal(p2, c(-1, 1, 1))
  

})




test_that('mirrorfill returns errors when it should.', {

  expect_error(mirrorfill(s, l1 = c(2:7, 9:14), l2 = c(99, 100)))

})


test_that('mirrorfill works correctly.', {

  # Create array in z = 0 plane
  a <- array(rep(1:36, 4), dim = c(12, 3, 4))
  a[, 3, ] <- 0

  # Make points 1:4 not on plane and symetric
  a[1:4, 3, ] <- c(1, -1, 2, -2)
  a[c(1, 3), 1:2, ] <- a[c(2, 4), 1:2, ]

  # Remove some data points
  missinga <- a
  missinga[c(1, 3), , 1] <- NA

  mirrorA <- mirrorfill(missinga, l1 = 5:12, l2 = 1:4)

  expect_equal(mirrorA, a)

  # Check l2 as a matrix works
  l2.m <- matrix(1:4, byrow = TRUE, ncol = 2)
  mirrorA2 <- mirrorfill(missinga, l1 = 5:12, l2 = l2.m)

  expect_equal(mirrorA2, a)


})



test_that('mirrorfill works correctly when multiple specimens have missing data.', {

  # Create array in z = 0 plane
  a <- array(rep(1:36, 4), dim = c(12, 3, 4))
  a[, 3, ] <- 0

  # Make points 1:4 not on plane and symetric
  a[1:4, 3, ] <- c(1, -1, 2, -2)
  a[c(1, 3), 1:2, ] <- a[c(2, 4), 1:2, ]

  # Remove some data points
  missinga <- a
  missinga[c(1, 3), , 1:3] <- NA


  mirrorA <- mirrorfill(missinga, l1 = 5:12, l2 = 1:4)

  expect_equal(mirrorA, a)

  # Check l2 as a matrix works
  l2.m <- matrix(1:4, byrow = TRUE, ncol = 2)
  mirrorA2 <- mirrorfill(missinga, l1 = 5:12, l2 = l2.m)

  expect_equal(mirrorA2, a)
})







test_that('mirror fill handles all missing vs not missing cases correctly.', {
    # In list of integers a, b, c, d,
    #   a and b are a pair and c and d are a pair
    #   Compare each pair both ways and replace missing values with mirrored version
    #   i.e. is a has missing data, replace a with reflection of b
    #   but if b has missing data, replace b with reflection of a
    #   If neither has missing data do nothing. If both have missing data, do nothing.


  # Create array in z = 0 plane
  a <- array(rep(1:36, 4), dim = c(12, 3, 4))
  a[, 3, ] <- 0

  # Make points 1:4 not on plane and symetric
  a[1:8, 3, ] <- c(1, -1, 2, -2)
  a[c(1, 3, 5, 7), 1:2, ] <- a[c(2, 4, 6, 8), 1:2, ]


  # Make a point that DOESN't reflect properly. 
  # This will be used to check that mirrorfill doesn't overwrite existing data.
  a[7, 1, 4] <- 9999

  # Remove some data points
  missinga <- a


  # 1 AND 2 missing.
  # 3, NOT 4 missing
  # (Not 5) AND 6 missing
  # NOT (4 OR 14) missing

  missinga[c(1, 2), , 1] <- NA
  missinga[3, , 2] <- NA
  missinga[6, , 3] <- NA

  
  mirrorA <- mirrorfill(missinga, l1 = 9:12, l2 = 1:8)

  # 1 AND 2 missing. Should still be NAs
  expect_true(all(is.na(mirrorA[c(1, 2), , 1])))

  # 3, NOT 4 missing
  # (Not 3) AND 5 missing
  # Should equal original a array
  expect_equal(mirrorA[, , 2], a[, , 2])
  expect_equal(mirrorA[, , 3], a[, , 3])

  # NOT 4 OR 14 missing
  #   mirrorA[7, 4, 1] should be 9999 NOT mirrorA[1a, 4, 1] 

  expect_equal(mirrorA[7, 1, 4], 9999)
  expect_false(mirrorA[7, 1, 4] ==  mirrorA[4, 1, 3])

  # Just to illustrate. 
  # Other mirrored points will be equal to their equivalent 
  #   point in other speciments
  expect_true(mirrorA[4, 1, 3] ==  mirrorA[4, 1, 2])
  expect_true(mirrorA[4, 1, 3] ==  mirrorA[4, 1, 1])
})
  



test_that('Mirrorfill example works', {
  
  
  # Make objects that are in the z = 0 plane.
  A <- array(rnorm(12 * 3 * 4), dim = c(12, 3, 4))
  
  # Make it symmetrical
  A[, 3, ] <- 0.1
  A2 <- A
  A2[, 3, ] <- -0.1
  A <- abind::abind(A, A2, along = 1)
  
  
  # Add some missing data points
  missinga <- A
  missinga <- abind::abind(missinga, array(NA, dim = c(2, 3, 4)), along = 1)
  
  mirrorA <- mirrorfill(missinga, l1 = c(1:24), l2 = c(23, 25, 24, 26))
  expect_equal(as.vector(A[23:24, 1:2,]),  as.vector(mirrorA[25:26, 1:2,]))
  expect_true(all(mirrorA[25:26, 3,] == 0.1))
})
