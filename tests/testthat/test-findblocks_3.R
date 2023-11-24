test_that("findblocks_testing1", {
  ## First test - No missing points
  mask <- rep(TRUE, 5)
  result <- findBlocks(mask)
  expect_equal(result, NULL)
  mask2 <- logical(0)
  expect_error(findBlocks(mask2), "Input 'mask' is empty")
  mask3 <- rep(TRUE, 10)
  result2 <- findBlocks(mask)
  expect_equal(result, NULL)
  mask4 <- c(1, 2, 3, 4, 5)
  expect_error(findBlocks(mask4), "Input 'mask4' must be a logical vector")
})
