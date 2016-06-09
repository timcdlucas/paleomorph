context('Test all functions used for mirroring')



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

  # Make an object that is reflected in the horizontal, z = 0 plane.
  s <- matrix(rep(1:21, 2), byrow = TRUE, ncol = 3)
  s[1:7, 3] <- -s[1:7, 3]


  pl <- bestplane(s)

  expect_equal(pl$n, c(0,0,1))  
  expect_equal(drop(pl$n[1] + pl$n[2] - pl$d), 0.0)


})

test_that('midline returns best plane', {

  # Make an object that is reflected in the horizontal, z = 0 plane.
  s <- matrix(rep(1:21, 2), byrow = TRUE, ncol = 3)
  s[1:7, 3] <- -s[1:7, 3]

  pl <- midline(s, 1:14)
  pl2 <- midline(s, c(2:5, 9:12))

  expect_equal(drop(pl$n[1] + pl$n[2] - pl$d), 0.0)
  expect_equal(drop(pl2$n[1] + pl2$n[2] - pl2$d), 0.0)

})



test_that('bestplane and midline return same values when no missing data.', {


  # Make an object that is reflected in the horizontal, z = 0 plane.
  s <- matrix(rep(1:21, 2), byrow = TRUE, ncol = 3)
  s[1:7, 3] <- -s[1:7, 3]

  p1 <- bestplane(s)
  l1 <- 1:14
  p2 <- midline(s, l1)

  expect_equal(p1, p2)



  set.seed(99)
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

  expect_equal(mirrorS, s)

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




test_that('mirrorfill returns errors when it should.', {

  expect_error(mirrorfill(s, l1 = c(2:7, 9:14), l2 = c(99, 100)))

})


test_that('mirrorfill works correctly.', {

  # Create array
  a <- array(rep(1:36, each = 4), dim = c(12, 3, 4))

  # Make it symmetric
  a[7:12, 1:2, ] <- a[1:6, 1:2, ]
  a[7:12, 3, ] <- -a[1:6, 3, ]

  # Remove some data points
  missinga <- a
  missinga[1:2, , 1] <- NA

  mirrorA <- mirrorfill(missinga, l1 = c(3:6, 9:12), l2 = c(1, 7, 2, 8))

  expect_equal(mirrorA, a)


})



test_that('mirrorfill works correctly when multiple specimens have missing data.', {

  # Create array
  a <- array(rep(1:36, by = 4), dim = c(12, 3, 4))

  # Make it symmetric
  a[7:12, 1:2, ] <- a[1:6, 1:2, ]
  a[7:12, 3, ] <- -a[1:6, 3, ]

  # Remove some data points
  missinga <- a
  missinga[1:2, , 1:3] <- NA

  mirrorA <- mirrorfill(missinga, l1 = c(3:6, 9:12), l2 = c(1, 7, 2, 8))

  expect_equal(mirrorA, a)

})







test_that('mirror fill handles all missing vs not missing cases correctly.', {
    # In list of integers a, b, c, d,
    #   a and b are a pair and c and d are a pair
    #   Compare each pair both ways and replace missing values with mirrored version
    #   i.e. is a has missing data, replace a with reflection of b
    #   but if b has missing data, replace b with reflection of a
    #   If neither has missing data do nothing. If both have missing data, do nothing.


  # Create array
  a <- array(rep(1:60), dim = c(20, 3, 4))

  # Make it symmetric
  a[11:20, 1:2, ] <- a[1:10, 1:2, ]
  a[11:20, 3, ] <- -a[1:10, 3, ]

  # Make a point that DOESN't reflect properly. 
  # This will be used to check that mirrorfill doesn't overwrite existing data.
  a[4, 1, 4] <- 9999

  # Remove some data points
  missinga <- a
  # 1 AND 11 missing.
  # 2, NOT 12 missing
  # Not 3 AND 13 missing
  # NOT 4 OR 14 missing

  missinga[c(1, 11), , 1] <- NA
  missinga[2, , 2] <- NA
  missinga[13, , 3] <- NA
  
  mirrorA <- mirrorfill(missinga, l1 = c(5:10, 15:20), l2 = c(1, 11, 2, 12, 3, 13, 3, 14))

  # 1 AND 11 missing. Should still be NAs
  expect_true(all(is.na(mirrorA[c(1, 11), , 1])))

  # 2, NOT 12 missing
  # Not 3 AND 13 missing
  # Should equal original a array
  expect_equal(mirrorA[, , 2], a[, , 2])
  expect_equal(mirrorA[, , 3], a[, , 3])

  # NOT 4 OR 14 missing
  #   mirrorA[4, 4, 1] should be 9999 NOT mirrorA[1a, 4, 1] 

  expect_equal(mirrorA[4, 1, 4], 9999)
  expect_false(mirrorA[4, 1, 4] ==  mirrorA[4, 1, 3])

  # Just to illustrate. 
  # Other mirrored points will be equal to their equivalent 
  #   point in other speciments
  expect_true(mirrorA[4, 1, 3] ==  mirrorA[4, 1, 2])
  expect_true(mirrorA[4, 1, 3] ==  mirrorA[4, 1, 1])
})


