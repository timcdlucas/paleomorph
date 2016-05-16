context('Test that EMMLi works correctly.')


test_that('EMMLi runs and creates output.', {

  set.seed(1)
  file <- paste0(tempdir(), 'EMMLiTest.csv')
  
  dat <- matrix(runif(36, -1, 1), ncol = 6, nrow = 6)
  diag(dat) <- 1

  mod1 <- data.frame(landmarks = letters[1:6], 
               modelA = rep(c(1, 2), each = 3),
               modelB = rep(c(1,2), 3),
               modelC = rep(c(1:3), 2)) 


  EMMLi(dat, 20, mod1, file)

  expect_true(file.exists(file))
  
  xx <- read.csv(file)

  expect_true(exists('xx'))

  # Bit of a mess because the csv is multiple tables.
  # So extract the top table only
  top <- xx[1:(which(xx[, 1] == '')[1] - 2), ]

  expect_true(sum(top$dAICc == 0) == 1)
  expect_true(all(as.character(top$dAICc) >= 0))
  expect_true(all(as.character(top$Model_l) >= 0))
  expect_true(all(as.character(top$Post_Pob) >= 0))
  expect_true(all(as.character(top$AICc) >= 0))

  unlink(file)

})


test_that('EMMLi can take mod with numbers in columns.', {

  set.seed(1)
  file <- paste0(tempdir(), 'EMMLiTest.csv')
  
  dat <- matrix(runif(36, -1, 1), ncol = 6, nrow = 6)
  diag(dat) <- 1

  mod1 <- data.frame(landmarks = letters[1:6], 
               model1 = rep(c(1, 2), each = 3),
               model2 = rep(c(1,2), 3),
               model3= rep(c(1:3), 2)) 



  EMMLi(dat, 20, mod1, file)

  expect_true(file.exists(file))
  
  xx <- read.csv(file)

  expect_true(exists('xx'))

  # Bit of a mess because the csv is multiple tables.
  # So extract the top table only
  top <- xx[1:(which(xx[, 1] == '')[1] - 2), ]

  expect_true(sum(top$dAICc == 0) == 1)


  unlink(file)


})




test_that('EMMLi is invariant to model order.', {

  set.seed(1)
  file <- paste0(tempdir(), 'EMMLiTest.csv')
  
  dat <- matrix(runif(36, -1, 1), ncol = 6, nrow = 6)
  diag(dat) <- 1

  mod1 <- data.frame(landmarks = letters[1:6], 
               model1 = rep(c(1, 2), each = 3),
               model2 = rep(c(2, 1), each = 3)) 



  EMMLi(dat, 20, mod1, file)

  expect_true(file.exists(file))
  
  xx <- read.csv(file)

  expect_true(exists('xx'))

  # Bit of a mess because the csv is multiple tables.
  # So extract the top table only
  top <- xx[1:(which(xx[, 1] == '')[1] - 2), ]

  expect_true(length(unique(top$MaxL)) == 3)
  expect_true(length(unique(top$dAICc)) == 3)

  unlink(file)


})





test_that('EMMLi can take mod with NAs.', {

  set.seed(1)
  file <- paste0(tempdir(), 'EMMLiTest.csv')
  
  dat <- matrix(runif(36, -1, 1), ncol = 6, nrow = 6)
  diag(dat) <- 1

  mod1 <- data.frame(landmarks = letters[1:6], 
               model1 = rep(c(1, 2), each = 3),
               model2 = rep(c(1,2), 3),
               model3= c(rep(c(1,2), 2), NA, NA)) 



  EMMLi(dat, 20, mod1, file)

  expect_true(file.exists(file))
  
  xx <- read.csv(file)

  expect_true(exists('xx'))

  # Bit of a mess because the csv is multiple tables.
  # So extract the top table only
  top <- xx[1:(which(xx[, 1] == '')[1] - 2), ]

  expect_true(sum(top$dAICc == 0) == 1)


  unlink(file)

})






