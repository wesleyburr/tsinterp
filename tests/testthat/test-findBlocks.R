test_that("Test for findBlocks function", {
  data("flux")
  z1 <- flux$SagOrig[1:100]
  z2 <- flux$PentOrig[1:100]
  z3 <- flux$S[1:100]
  z1[which(z3 == FALSE)] <- NA
  blocks <- findBlocks(z1)
  expect_equal(class(blocks), "numeric")
  expect_equal(blocks[1], 89)
  expect_equal(findBlocks(z2), NULL)
  expect_length(blocks, 3)
})
