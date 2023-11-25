## Need to add code to findPowers for error (stopifnot), like putting in negative n value.

test_that("findPowers_test", { 
  expect_error(findPower(-1,1,1,1), "Error")
})
