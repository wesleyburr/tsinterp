# test 1 for linInt() 
test_that("test to see if linInt() work",{
  set.seed(2)
  x <- runif(20,1,10)
  miss_block <- c(2,3,7,8,9,18,19)
  x[miss_block] <- NA
  blocks <- findBlocks(x)
  result <- c( 2.66394034, 2.61344932, 2.56295830, 2.51246728, 9.49455405, 9.49127463, 
              8.60591939, 7.72056415, 6.83520891, 5.94985368, 5.97406660, 3.15005284, 
              7.84461982, 2.62738091, 4.64753963, 8.68193608, 9.787586, 7.083329, 
              4.379072, 1.674815)
  expect_equal(linInt(x,blocks),result)
})


# test 2 for linInt(), aiming at error detection
test_that("test for error triggering in linInt",{
  expect_error(linInt(NULL))
  expect_error(linInt("abc"))
  expect_error(linInt(TRUE))
  expect_error(linInt(data.frame(runif(10))))
})
  