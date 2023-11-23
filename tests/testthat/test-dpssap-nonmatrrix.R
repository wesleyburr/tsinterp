test_that("dpssap_test", { 
  m <- 1
  expect_error(dpssap(m, 1), "is.matrix(V) is not TRUE")
})
