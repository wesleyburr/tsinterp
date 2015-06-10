summaryRprof(filename = "1-estimateMT.out")


runExample <- function() {
  data("flux")
  
  z1 <- flux$SagOrig
  z1[which(flux$S == FALSE)] <- NA
  
  # Unfortunately, not fast enough to run for CRAN checks
  sagInt <- interpolate(z = z1, gap = which(flux$S == FALSE), maxit = 3, delT = 86400)
  
}


Rprof("3-interpolate.out")
runExample()
Rprof(NULL)

Rprof("2-bivarinter.out")
z1 <- flux$SagOrig
z1[which(flux$S == FALSE)] <- NA
z2 <- flux$PentOrig
 sagInt <- BiVarInt(z1 = z1, z2 = z2, gap1 = which(flux$S == FALSE), 
                   gap2 = NULL, maxit = 3, delT = 86400)

Rprof(NULL)


############# Compare results to before changing code 

data("flux")

z1 <- flux$SagOrig
z1[which(flux$S == FALSE)] <- NA

# Unfortunately, not fast enough to run for CRAN checks
sagInt2 <- interpolate(z = z1, gap = which(flux$S == FALSE), maxit = 3, delT = 86400)

originalIntZf <- scan(file = "originalInter/originalZf-interpolate.txt", sep = ";")
originalIntP <- scan(file = "originalInter/originalP-interpolate.txt")
originalIntConverged  <- scan(file = "originalInter/originalConverged-interpolate.txt", what = logical())
originalIntDiffC  <- scan(file = "originalInter/originalDiffC-interpolate.txt")


all.equal(sagInt[[1]], originalIntZf);

