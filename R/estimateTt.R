#' Estimate the sinusoidal components in a time series.
#'
#' This function estimates sinusoidal components in a time series using the multitaper method.
#'
#' @param x The input time series.
#' @param epsilon A tolerance for peak estimation convergence.
#' @param dT The time interval between data points.
#' @param nw The time-half bandwidth parameter for multitaper.
#' @param k The number of tapers.
#' @param sigClip Significance level for peak detection.
#' @param progress If TRUE, show progress information.
#' @param freqIn A pre-defined set of frequencies (optional).
#'
#' @return A matrix containing sinusoidal components.
#'
#' @details
#' The function follows these algorithmic steps:
#' 1. Spectrum/F-test pilot estimate using the multitaper method.
#' 2. Estimation of significant peaks based on the significance level (sigClip).
#' 3. Refinement of peak locations by eliminating duplicates and optimizing peak frequencies.
#' 4. Estimation of phase and amplitude by inverting the spectrum (line component removal).
#'
#' @export
#' 
#' @importFrom multitaper spec.mtm
#' @examples
#' estimateTt(x, epsilon = 1e-10, dT = 0.01, nw = 4, k = 5, sigClip = 0.05)
#'
estimateTt <- function(x, epsilon, dT, nw, k, sigClip, progress=FALSE, freqIn=NULL) {
  
  ################################################################################
  # Algorithm step 1: spectrum/Ftest pilot estimate
  pilot <- multitaper::spec.mtm(x, deltat=dT, nw=nw, k=k, Ftest=TRUE, plot=FALSE)
  
  if(is.null(freqIn)) {
    ################################################################################
    # Algorithm step 2: estimate significant peaks (sigClip)
    fmesh <- pilot$mtm$Ftest
    fsig <- fmesh > qf(sigClip, 2, pilot$mtm$k)
    floc <- which(fsig==TRUE)
    if(length(floc > 0)) {
      ###########################################################################   
      delta <- floc[2:length(floc)] - floc[1:(length(floc)-1)]
      if(length(which(delta==1)) > 0) {
        bad <- which(delta==1)
        if(!is.null(bad)) {
          if(progress) {
            for(j in 1:length(bad)) {
              cat(paste("Peak at ", formatC(pilot$freq[floc[bad[j]]], width=6, format="f"),
                        "Hz is smeared across more than 1 bin. \n", sep=""))
            } 
          }
        }
        floc <- floc[-bad] # eliminate the duplicates
      }
      
      ################################################################################
      # Algorithm step 3: estimate centers
      dFI <- pilot$freq[2]
      # epsilon <- 1e-10
      maxFFT <- 1e20
      max2 <- log(maxFFT, base=2)
      max3 <- log(maxFFT, base=3)
      max5 <- log(maxFFT, base=5)
      max7 <- log(maxFFT, base=7)
      
      freqFinal <- matrix(data=0, nrow=length(floc), ncol=1)
      
      for(j in 1:length(floc)) {
        if(progress) {
          cat(".")  
        }
        f0 <- pilot$freq[floc[j]]
        
        if(progress) {
          cat(paste("Optimizing Peak Near Frequency ", f0, "\n", sep=""))
        }
        
        # increasing powers of 2,3,5,7 on nFFT until peak estimate converges
        pwrOrig <- floor(log2(pilot$mtm$nFFT)) + 1
        fI <- f0
        converge <- FALSE
        
        for(k7 in 0:max7) {
          for(k5 in 0:max5) {
            for(k3 in 0:max3) {
              for(k2 in 1:max2) {
                
                if(!converge) {
                  nFFT <- 2^pwrOrig * 2^k2 * 3^k3 * 5^k5 * 7^k7
                  tmpSpec <- multitaper::spec.mtm(x, deltat=dT, nw=5, k=8, plot=FALSE, Ftest=TRUE,
                                      nFFT=nFFT)
                  dF <- tmpSpec$freq[2]
                  f0loc <- which(abs(tmpSpec$freq - f0) <= dF)
                  range <- which(tmpSpec$freq <= (f0+1.1*dFI) & tmpSpec$freq >= (f0-1.1*dFI))
                  
                  fI2 <- tmpSpec$freq[which(tmpSpec$mtm$Ftest == max(tmpSpec$mtm$Ftest[range]))]
                  if(abs(fI - fI2) > epsilon) {
                    fI <- fI2
                  } else {
                    fF <- fI2
                    converge <- TRUE
                  }
                }
              }}}}
        freqFinal[j] <- fF
        if(progress) {
          cat(paste("Final frequency estimate: ", fF, "\n", sep=""))
        }
      }
      if(progress) {
        cat("\n")
      }
    } else {
      freqFinal <- NULL
      floc <- -1
    } # end of "there are freqs detected"
  } else {  # case where frequencies are already obtained
    freqFinal <- freqIn
    floc <- 1:length(freqFinal)
    if(length(freqFinal)==1 & freqFinal[1]==0) {
      floc <- -1
    }
  }
  ################################################################################
  # Algorithm step 4: frequencies obtained, estimate phase and amplitude
  #    by inverting the spectrum (i.e. line component removal)
  if(length(floc) > 1 | floc[1] > 0) {
    sinusoids <- matrix(data=0, nrow=length(x), ncol=length(floc))
    amp <- matrix(data=0, nrow=length(floc), ncol=1)
    phse <- matrix(data=0, nrow=length(floc), ncol=1)
    N <- length(x)
    t <- seq(1, N*dT, dT)
    
    for(j in 1:length(floc)) {
      sinusoids[, j] <- removePeriod(x, freqFinal[j], nw=5, k=8, deltaT=dT, warn=FALSE, prec=1e-10, sigClip=sigClip) 
      fit <- lm(sinusoids[, j] ~ sin(2*pi*freqFinal[j]*t) + cos(2*pi*freqFinal[j]*t) - 1)
      phse[j] <- atan(fit$coef[2] / fit$coef[1])
      amp[j] <- fit$coef[1] / cos(phse[j])
    }
    
    attr(sinusoids, "Phase") <- phse
    attr(sinusoids, "Amplitude") <- amp
    attr(sinusoids, "Frequency") <- freqFinal
    return(sinusoids)
  } else {
    sinusoids <- matrix(data=0, nrow=length(x), ncol=length(floc))
    attr(sinusoids, "Phase") <- 0
    attr(sinusoids, "Amplitude") <- 0
    attr(sinusoids, "Frequency") <- 0
    return(sinusoids)
  }
}
