#' Utility that takes pre-interpolated (linear?) series and computes spectra
#'and ACVF / CCVF series from them, which are then passed to mwXSwiener.
#'
#' @param xd1 
#' @param xd2 
#' @param ok1 
#' @param ok2 
#' @param dT 
#' @param blocks 
#' @param neh 
#' @param maxlag 
#' @param clipMax 
#'
#' @return
#' @export
#'
#' @examples
"estimateWt" <- function(xd1, xd2, ok1, ok2, dT, blocks, neh, maxlag, clipMax) {
  
  if(length(which(is.na(xd1))) > 0 | length(which(is.na(xd2))) > 0) {
    stop( "Series must be pre-interpolated before passing to WtEstimator." )
  }
  # need the spectrum to be big enough to provide all the lags necessary ...
  nFFT <- 2^(floor(log2(3*maxlag))+2)
  
  sp1 <- spec.mtm(xd1, deltat=dT, nFFT=nFFT, plot=FALSE, returnInternals=TRUE)
  sp2 <- spec.mtm(xd2, deltat=dT, nFFT=nFFT, plot=FALSE, returnInternals=TRUE)
  R11 <- SpecToACVdual(sp1, maxlag=2*maxlag)
  R22 <- SpecToACVdual(sp2, maxlag=2*maxlag)
  R12 <- spec.mtm.cross(sp1, sp2, maxlag=2*maxlag)
  R21 <- spec.mtm.cross(sp2, sp1, maxlag=2*maxlag)
  
  # loop on gaps in xd1
  yd1 <- xd1
  # cat("starting \n")
  
  
  # ************** ERROR HERE
  #  ** if neh > blocks[m, 1] this will crash
  #  ** wants lots of 'good' data to its left
  
  for(m in 1:length(blocks[, 1]))  {
    rng <- max(1, (blocks[m, 1] - neh)):(blocks[m, 2] + neh)
    fill <- mwXSwiener(xd1[rng], xd2[rng], ok1[rng], ok2[rng], R11, R12, R21, R22)
    fill[which(abs(fill) > clipMax)] <- sign(fill[which(abs(fill) > clipMax)])*clipMax
    yd1[blocks[m, 1]:blocks[m, 2]] <- fill[(neh+1):(neh+blocks[m, 3])]
  }
  # cat("\n")
  yd1
}
