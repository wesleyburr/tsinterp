#' Compute cross-covariance sequence from two spectral objects.
#'
#' This function computes the cross-covariance sequence from two spectral objects.
#' It requires that the spectral objects were computed with the `returnInternals` argument set to TRUE.
#' The two spectral objects must have the same parameters for successful computation.
#' The function can also limit the maximum lag for the cross-covariance sequence.
#'
#' @param sp1 The first spectral object.
#' @param sp2 The second spectral object.
#' @param maxlag Maximum lag for the cross-covariance sequence (optional).
# 
#' @return A cross-covariance sequence.
#'
#' @export
#' 
#' @importFrom multitaper spec.mtm
#' 
#' @details
#' The function computes the cross-covariance sequence by utilizing the spectral properties
#' of two input spectral objects. The `maxlag` parameter controls the maximum lag
#' for the cross-covariance sequence. Ensure that the spectral objects are computed with
#' the `returnInternals` argument set to TRUE and have consistent parameters.
# 
#' @examples
#' sp1 <- multitaper::spec.mtm(x1, nw = 4, k = 5, returnInternals = TRUE)
#' sp2 <- multitaper::spec.mtm(x2, nw = 4, k = 5, returnInternals = TRUE)
#' spec.mtm.cross(sp1, sp2)
#'
"spec.mtm.cross" <- function(sp1, sp2, maxlag=NULL) {
  
  # make sure objects were returned "withInternals"
  if(is.null(sp1$mtm$eigenCoefs) | is.null(sp2$mtm$eigenCoefs)) { stop("Spectrum objects must have been computed with returnInternals=TRUE.") }
  if(sp1$freq[2] != sp2$freq[2]) { stop("dF not the same for two series.") }
  if(sp1$mtm$nfreqs != sp2$mtm$nfreqs) { stop("nFFT not the same for the two series.") }
  if(sp1$mtm$k != sp2$mtm$k) { stop("Number of tapers not the same for the two series.") }
  if(is.null(maxlag)) { maxlag = sp1$mtm$nfreqs }
  
  dF = sp1$freq[2]
  cFreq <- sp1$mtm$nfreqs
  K <- sp1$mtm$k
  
  X8 <- 1/K * (rowSums(Re(sp1$mtm$eigenCoefs) * Re(sp2$mtm$eigenCoefs)) + 
                 rowSums(Im(sp1$mtm$eigenCoefs) * Im(sp2$mtm$eigenCoefs)) )
  
  Y8 <- 1/K * (rowSums(Im(sp1$mtm$eigenCoefs) * Re(sp2$mtm$eigenCoefs)) -
                 rowSums(Re(sp1$mtm$eigenCoefs) * Im(sp2$mtm$eigenCoefs)) )
  s12 <- complex(real=X8, imaginary=Y8)
  X8 <- X8 * dF
  Y8 <- Y8 * dF
  X8 <- c(X8, rev(X8[-c(1, cFreq)]))
  Y8 <- c(Y8, -1*rev(Y8[-c(1, cFreq)]))
  r12 <- Re(fft(complex(real=X8,imaginary=Y8)))
  r12 <- c(r12[(cFreq+1):length(r12)], r12[1:cFreq])[(cFreq-1-maxlag):(cFreq-1+maxlag)]
  attr(r12, "ZeroF") <- maxlag+1
  attr(r12, "MaxLag") <- maxlag
  return(r12)
}
