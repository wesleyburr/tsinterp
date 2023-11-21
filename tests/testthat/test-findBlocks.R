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

