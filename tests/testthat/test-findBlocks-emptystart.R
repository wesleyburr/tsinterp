## This test fails... issue detecting missing block when it is at the start

test_that("Testing findBlocks",{
  data("flux")
  missing_mnt <- flux$Mnt
  missing_mnt[missing_mnt == 5] <- NA
  blocks <- findBlocks(missing_mnt)
  expect_equal(blocks[1,1], 1)
})
