test_that("findBlocks_test", {
  data("flux")
  miss <- flux$S
  miss[miss == FALSE] <- NA
  blocks <- findBlocks(miss)
  expect_equal(ncol(blocks), 3)
  expect_equal(class(blocks)[1], "matrix")
  expect_equal(class(blocks)[2], "array")
  expect_length(blocks, 1422)
})
