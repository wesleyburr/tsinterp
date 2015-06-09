summaryRprof(filename = "1-estimateMT.out")


runExample <- function() {
  data("flux")
  
  z1 <- flux$SagOrig
  z1[which(flux$S == FALSE)] <- NA
  
  # Unfortunately, not fast enough to run for CRAN checks
  sagInt <- interpolate(z = z1, gap = which(flux$S == FALSE), maxit = 3, delT = 86400)
  
}


Rprof("2-interpolate.out")
runExample()
Rprof(NULL)



