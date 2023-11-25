test_that("Test for linear Interpolation function", {
  data(flux)
  miss <- flux$S
  miss[miss == FALSE] <- NA
  blocks <- findBlocks(miss)
  fluxInt <- linInt(flux$SagOrig, blocks)
  expect_equal(class(fluxInt), "numeric")
  expect_error(linInt(flux))
  expect_equal(length(fluxInt), length(flux$SagOrig))
})
