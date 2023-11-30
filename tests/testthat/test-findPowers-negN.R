## Fails, check error message 

test_that("findPowers_test", { 
  expect_error(findPower(-1,1,1,1), "N must be greater than or equal to 0")
})
