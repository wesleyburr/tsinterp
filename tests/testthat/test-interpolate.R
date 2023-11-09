test_that("interpolate flux", {
  source("C:/Users/Pars.R/Documents/tsinterp/R/Interpolate.R")
  load("flux.RData")
  
  z1 <- flux$SagOrig
  z1[which(flux$S == FALSE)] <- NA
  interpolate(z = z1, gap = which(flux$S == FALSE), maxit = 3, delT = 86400)
  expect_equal(interpolate(z = z1, gap = which(flux$S == FALSE),
                           maxit = 3, delT = 86400)
               , load("sagInt.RData"))
})
