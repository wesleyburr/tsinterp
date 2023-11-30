test_that("test for findBlocks works", {
  set.seed(2)
  x <- runif(20,1,10)
  miss_block <- c(2,3,7,8,9,19,20)
  # result1 has the block at the end of the TS
  # result2 doesn't have the block at the end of the TS
  result1 <- matrix(c(2,3,2,7,9,3,19,20,2),byrow=TRUE,ncol=3)
  result2 <- matrix(c(2,3,2,7,9,3),byrow=TRUE,ncol =3)
  x[miss_block] <- NA
  expect_type(findBlocks(x),'matrix')
  expect_equal(findBlocks(x), result2)
  x <- runif(20,1,10)
  expect_equal(findBlocks(x), NULL)
})
# second test for findBlocks(), aiming at error detection
test_that("test for error triggering in findBlocks",{
  expect_error(findBlocks(NULL))
  expect_error(findBlocks("abc"))
  expect_error(findBlocks(TRUE))
  expect_error(findBlocks(data.frame(runif(10))))
})

# test 3 for findBlocks() about warnings
test_that("test for start warning message in findBlocks",{
  set.seed(2)
  x <- runif(20,1,10)
  miss_block <- c(1,2,3,7,8,9,19)
  x[miss_block] <- NA
  # time series x has missing the block at the start
  expect_warning(findBlocks(x))
})

# test 4 for findBlocks() about warnings
test_that("test for end warning message in findBlocks",{
  set.seed(2)
  x <- runif(20,1,10)
  miss_block <- c(2,3,7,8,9,19,20)
  x[miss_block] <- NA
  # time series x has missing the block at the end
  expect_warning(findBlocks(x))
})

# test 5 for findBlocks() about no error
test_that("test for no error message in findBlocks",{
  set.seed(2)
  x <- runif(20,1,10)
  miss_block <- c(2,3,7,8,9,18,19)
  x[miss_block] <- NA
  # time series x has the right form for findBlocks()
  expect_no_error(findBlocks(x))
})

# test 6 for findBlocks() about no warnings
test_that("test no warning message in findBlocks",{
  set.seed(2)
  x <- runif(20,1,10)
  miss_block <- c(2,3,7,8,9,18,19)
  x[miss_block] <- NA
  # time series x has the right form for findBlocks()
  expect_no_warning(findBlocks(x))
})




