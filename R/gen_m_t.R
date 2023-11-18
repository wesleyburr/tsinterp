#' @title gen_m_t is a used to generate mt for the interpolation function
#' 
#' @details A function that, given argument zI, computes a MtP/TtP
#' pairing of estimates, which can then be used in the 
#' maximum-likelihood iterative step of the EM algorithm
#' for interpolation.
#'
#' @param zI   Input time series; should have no missing values.
#' @param N    Length of zI
#' @param delT   delta-T; the time step of the time series.
#' @param sigClip   Significance level for periodic component detection.
#' @param progress  Should a progress tick be printed or not?
#'
#' @return a class list that has Ttp and freqsave as objects
#' 
#' 
#' 
#' 
#' @export
#'
#' 
#' 
gen_m_t <- function(zI, N, delT, sigClip, progress) {
  MtP <- estimateMt(x=zI, N=N, nw=5, k=8, pMax=2)
  TtTmp <- estimateTt(x=zI - MtP, epsilon=1e-6, dT=delT, nw=5, k=8,
                      sigClip=0.05, progress=FALSE)
  freqRet <- attr(TtTmp, "Frequency")
  if(length(freqRet) > 1 | (length(freqRet)==1 && freqRet != 0)) {
    TtP <- rowSums(TtTmp) 
  } else {
    TtP <- TtTmp
  }
  freqSave <- attr(TtTmp, "Frequency")
  # returns
  return(list(TtP = TtP,
              freqSave = freqSave
  )
  )
} 

