test_that("findblocks_testing1", {
  ## First test - No missing points
  mask <- rep(TRUE, 5)
  result <- findBlocks(mask)
  expect_equal(result, NULL)
})
