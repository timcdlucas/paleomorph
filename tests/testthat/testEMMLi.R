]context('Test that EMMLi works correctly.')


test_that('EMMLi runs and creates output.', {

  set.seed(1)
  file <- paste0(tempdir(), 'EMMLiTest.csv')
  
  dat <- matrix(runif(36, -1, 1), ncol = 6, nrow = 6)
  diag(dat) <- 1

  mod1 <- data.frame(landmarks = letters[1:6], 
               modelA = rep(c(1, 2), each = 3),
               modelB = rep(c(1,2), 3),
               modelC = rep(c(1:3), 2)) 

  #varList <- c('mod$modelA', 'mod$modelB', 'mod$modelC')


  EMMLi(dat, 20, mod1, file)

  expect_true(file.exists(file))
  
  xx <- read.csv(file)

  expect_true(exists('xx'))

  # Bit of a mess because the csv is multiple tables.
  # So extract the top table only
  top <- xx[1:(which(xx$dAICc == 'unintegrated') - 2), ]

  expect_true(sum(top$dAICc == 0) == 1)


  unlink(file)

})


test_that('EMMLi can take mod with numbers in columns.', {

  set.seed(2)
  file <- paste0(tempdir(), 'EMMLiTest.csv')
  
  dat <- matrix(runif(36, -1, 1), ncol = 6, nrow = 6)
  diag(dat) <- 1

  mod1 <- data.frame(landmarks = letters[1:6], 
               model1 = rep(c(1, 2), each = 3),
               model2 = rep(c(1,2), 3),
               model3= rep(c(1:3), 2)) 

  #varList <- c('mod$model1', 'mod$model2', 'mod$model3')


  EMMLi(dat, 20, mod1, file)

  expect_true(file.exists(file))
  
  xx <- read.csv(file)

  expect_true(exists('xx'))

  # Bit of a mess because the csv is multiple tables.
  # So extract the top table only
  top <- xx[1:(which(xx$dAICc == 'unintegrated') - 2), ]

  expect_true(sum(top$dAICc == 0) == 1)


  unlink(file)


})
