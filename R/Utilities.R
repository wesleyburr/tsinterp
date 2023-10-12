#######################################################################
#
#    WtEstimator
#
#    Wrapper function for mwXSwiener
#
#    Takes pre-interpolated (linear?) series and computes spectra
#    and ACVF / CCVF series from them, which are then passed to 
#    mwXSwiener.
#
#######################################################################
estimateMt <- function(x, N, nw, k, pMax) {
    V <- dpss(n=N, nw=5, k=8)$v
    test <- dpssap(V,pMax) # fit quadratic
    U <- test[[1]]
    R <- test[[2]]
    Y <- t(V) %*% x
    a <- t(U) %*% Y
    r <- Y - U %*% a
    xhat <- V %*% r + R %*% a
    phat <- R %*% a
    return(phat)
}

estimateTt <- function(x, epsilon, dT, nw, k, sigClip, progress=FALSE, freqIn=NULL) {

    ################################################################################
    # Algorithm step 1: spectrum/Ftest pilot estimate
    pilot <- spec.mtm(x, deltat=dT, nw=nw, k=k, Ftest=TRUE, plot=FALSE)

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
            tmpSpec <- spec.mtm(x, deltat=dT, nw=5, k=8, plot=FALSE, Ftest=TRUE,
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

findPowers <- function(N,f0,Nyq,prec) {
  nFFT <- 1e30

  low2 <- 0
  high2 <- floor(log(N,2))+2
  low3 <- 0
  high3 <- floor(log(N,3))+2
  low5 <- 0
  high5 <- floor(log(N,5))+2
  low7 <- 0
  high7 <- floor(log(N,7))+2
  for(i in low2:high2) {
    for(j in low3:high3) {
      for(k in low5:high5) {
        for(l in low7:high7) {
          att <- 2^i * 3^j * 5^k * 7^l
          if((att > 2*N) & att < 100*N) {
            df <- (Nyq*2)/att
            if( abs(trunc(f0/df)*df - f0) < prec ) {
              if(att < nFFT) {
                nFFT <- att 
              }
            }
          } # big enough
        } # end of 7
      } # end of 5
    } # end of 3
  } # end of 2
  if(nFFT == 1e30) {
    return(-1)
  } else {
    return(nFFT)
  }
}


#######################################################################
#
#   dpssap
#      
#   Computes the Discrete Prolate Spheroidal Sequences Associated 
#   Polynomials for given N and given (pre-computed) dw of dimension
#   N * K, concentration NW. Takes parameter maxdeg as maximum degree.
#
#   Adapted from 'dpssap.f', written by David J. Thomson, February 2001
#
#######################################################################
#   
#  Given:
#  x      Time indexes for evaluation of polynomials
#  V      Matrix of Slepian sequences, k=0,1,...,K-1, n=0,1,...,N-1
#  maxdeg Maximum Degree
#
#  Returns:
#  A      Orthonormal Projection Matrix, A[k,j] = Inner Product of
#           dw[t,k] with poly[t,j] on [1,N]. The polynomials are
#           defined by the constraint that
#             \sum_{k=0}^{K-1} A(k,i) * A(k,j) = \delta_{i,j}
#           Note that the polynomials are NOT orthogonal in the sense
#             \sum_{t=1}^{N} poly(t,i) * poly(t,j) \ne  \delta_{j,L}
#

dpssap <- function(V, maxdeg) {

    # Sanity checks
    stopifnot(is.matrix(V), is.numeric(maxdeg), maxdeg>=0)
    N <- length(V[, 1])
    K <- length(V[1, ])
    P <- maxdeg + 1
    timeArr <- 1:N

    R <- matrix(data=0, nrow=N, ncol=P)
    U <- matrix(data=0, nrow=K, ncol=P)

    # Setup centered time index
    midTime <- (1+N) / 2
    scl <- 2/(N-1)
    timeArrC <- (timeArr - midTime) * scl
    
    # Start with Gegenbauer polynomials; convergence is faster
    alpha <- 0.75
    R[, 1] <- 1.0
    if(maxdeg > 0) {
      R[, 2] <- 2 * alpha * timeArrC
      if(maxdeg > 1) {
        for(j in 2:maxdeg) {
          A1 <- 2 * ( (j-1) + alpha ) / j
          A2 <- ( (j-2) + 2 * alpha ) / j

          R[, (j+1)] <- A1 * timeArrC * R[, j] - A2 * R[, (j-1)]
        } # end of loop on higher orders
      } # end of maxdeg > 1
    } # end of maxdeg > 0

    # Inner Products of R and V
    for(L in 1:P) {
      Kmin <- ( (L-1) %% 2 ) + 1
      for(k in seq(Kmin, K, 2)) {  # loop on non-zero Slepians
        U[k, L] <- t(V[, k]) %*% R[, L]
      }
    }

    # Degree 0, 1 (manual) -- L = degree+1
    for(L in 1:min(2,P)) {
      scl <- 1 / sqrt( sum(U[, L]^2) )
      U[, L] <- U[, L] * scl # orthonormalize
      R[, L] <- R[, L] * scl
    }

    # loop on higher degrees, applying Gram-Schmidt only on similar
    # parity functions (as even/odd are already orthogonal in U)
    if( P > 2 ) {
      for(L in 3:P) {
        if(L %% 2 == 0) {
          Kmin <- 2
        } else {
          Kmin <- 1
        }
        for(j in seq(Kmin, L-1, 2)) {
          scl <- sum( U[, L] * U[, j] )
          U[, L] <- U[, L] - scl * U[, j] # Gram-Schmidt
          R[, L] <- R[, L] - scl * R[, j]
        }
        scl <- 1 / sqrt(sum(U[, L]^2))
        U[, L] <- U[, L] * scl  # orthonormalize
        R[, L] <- R[, L] * scl
      }
    }

    Hn <- colSums(R^2)
    return(list(U,R,Hn))
}

################################################################################
#
#  Find blocks of missing points
#
#  Input: mask of TRUE/NA (or anything/NA -- missing = NA only)
#
################################################################################
"findBlocks" <- function(mask) {
  Nlen <- length(mask)
  mask <- which(is.na(mask))
  # case: there are missing points
  if(length(mask) > 0) {
    diffs <- mask[-1] - mask[-length(mask)]
    diffs <- which(diffs > 1)

    # case: 1 gap only, possibly no gaps
    if(length(diffs)==0) {
      blocks <- matrix(data=0, nrow=1, ncol=3)
      blocks[1, 1:2] <- c(mask[1], mask[length(mask)])
    } else {
      blocks <- matrix(data=0,nrow=length(mask),ncol=3)
      blocks[1, 1] <- mask[1]
      blocks[1, 2] <- mask[diffs[1]]
      k <- 1
      for(j in 1:length(diffs)) {
        k <- k+1
        blocks[k, 1:2] <- c(mask[diffs[j]+1],mask[diffs[j+1]])
      }
      blocks[k,2] <- max(mask)
      blocks <- blocks[1:k, ]
    }
    blocks[,3] <- blocks[,2] - blocks[,1] + 1

    # checks to remove start/end of sequence
    if(blocks[1,1]==1) {
      blocks <- blocks[-1, ]
    }
    if(blocks[length(blocks[,1]),2]==Nlen) {
      blocks <- blocks[-length(blocks[,1]), ] 
    }
  } else {
    blocks <- NULL
  }
  blocks
}

######################################################################
#
#  Linear interpolator, requires data and blocks of missing points
#
######################################################################
"linInt" <- function(dat,blocks) {
  nGap <- length(blocks[,1])
  for(j in 1:nGap) {
    dY <- (dat[blocks[j,2]+1] - dat[blocks[j,1]-1])/(blocks[j,3]+1)
    st <- dat[blocks[j,1]-1]
    lt <- dat[blocks[j,2]+1]
    if(dY != 0) {
      fill <- seq(st,lt,dY)
    } else {
      fill <- rep(st,blocks[j,2]-blocks[j,1]+3) 
    }
    fill <- fill[c(-1,-length(fill))]
    dat[blocks[j,1]:blocks[j,2]] <- fill
  }
  dat
}

"SpecToACV" <- function(spec,maxlag) {
  s <- spec$spec
  dF <- spec$freq[2] 
  x <- matrix(data=0,nrow=(spec$mtm$nfreqs-1)*2,ncol=1)
  x[1:spec$mtm$nfreqs] = s*dF
  x[(spec$mtm$nfreqs+1):length(x)] <- x[(spec$mtm$nfreqs-1):2]
  x <- as.complex(x)
  x <- Re(fft(x,inverse=TRUE))
  x[1:(maxlag+1)]
}

#######################################################################
#
#    spec.mtm.cross
#
#    Computes cross-covariance sequence from two spectral objects
#
#######################################################################
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

#######################################################################
#
#    SpecToACV
#
#    Converts spectral object to ACVF sequence.
#
#######################################################################
"SpecToACVdual" <- function(spec,maxlag=NULL) {

  if(is.null(maxlag)) { maxlag = spec$mtm$nfreqs -2 }

  s <- spec$spec
  dF <- spec$freq[2] 
  x <- matrix(data=0,nrow=(spec$mtm$nfreqs-1)*2,ncol=1)
  x[1:spec$mtm$nfreqs] = s*dF
  x[(spec$mtm$nfreqs+1):length(x)] <- x[(spec$mtm$nfreqs-1):2]
  x <- as.complex(x)
  x <- Re(fft(x,inverse=TRUE))
  cFreq <- spec$mtm$nfreqs
  x <- c(x[(cFreq+1):length(x)], x[1:cFreq])[(cFreq-1-maxlag):(cFreq-1+maxlag)]
  attr(x, "ZeroF") <- maxlag+1
  attr(x, "MaxLag") <- maxlag
  x
}

