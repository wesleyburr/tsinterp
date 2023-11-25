## Test gets error, how to write properly?

test_that("dpssap_test", { 
  m <- matrix(1:9, nrow = 3, ncol = 3)
  expect_error(dpssap(m, -1), Error)
})
