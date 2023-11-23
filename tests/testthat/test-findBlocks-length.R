test_that("Testing findBlocks",{
  data("flux")
  missing_mnt <- flux$Mnt
  missing_mnt[missing_mnt == 5] <- NA
  blocks <- findBlocks(missing_mnt)
  expect_equal(blocks[,3], rep(31, nrow(blocks)))
  expect_equal(blocks[1,2]-blocks[1,1]+1, blocks[1,3])
})
