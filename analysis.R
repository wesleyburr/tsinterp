
runInterExample <- function(...) {
  data("flux")
  z1 <- flux$SagOrig
  z1[which(flux$S == FALSE)] <- NA
  # Unfortunately, not fast enough to run for CRAN checks
  sagInt <- interpolate(z = z1, gap = which(flux$S == FALSE), maxit = 3, delT = 86400, ...)
}

runBiVarExample <- function(...) {
  data("flux")
  
  z1 <- flux$SagOrig
  z1[which(flux$S == FALSE)] <- NA
  z2 <- flux$PentOrig
  
  # Unfortunately, not fast enough to run for CRAN checks
   sagInt <- BiVarInt(z1 = z1, z2 = z2, gap1 = which(flux$S == FALSE), 
                      gap2 = NULL, maxit = 3, delT = 86400, ...)
  
}

############# Compare results to before changing code 
# Unfortunately, not fast enough to run for CRAN checks
checkInterpolation <- function(...) {
  data("flux")
  z1 <- flux$SagOrig
  z1[which(flux$S == FALSE)] <- NA
  
  sagInt <- interpolate(z = z1, gap = which(flux$S == FALSE), maxit = 3, delT = 86400, ...)

  originalIntZf <- read.csv("../originalInter/originalZf-interpolate.csv")[ ,1]
  all.equal(sagInt[[1]], originalIntZf);
  print(ifelse(all.equal(sagInt[[1]], originalIntZf), "Result is ok", "result is wrong"))
  sagInt
}

######## Bi var Test
checkBiVar <- function(...) {
  data("flux")
  z1 <- flux$SagOrig
  z1[which(flux$S == FALSE)] <- NA
  z2 <- flux$PentOrig
  sagInt <- BiVarInt(z1 = z1, z2 = z2, gap1 = which(flux$S == FALSE), 
                    gap2 = NULL, maxit = 3, delT = 86400, ...)

  originalBivarZf <- read.csv("../originalBivar/originalZf-bivar.csv")[ , 1]
  print(ifelse(all.equal(sagInt[[1]], originalBivarZf), "Result is ok", "result is wrong"))
  
  sagInt
}