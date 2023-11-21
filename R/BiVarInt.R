#' @title BiVarInt
#' @param z1 first time series (possibly) with gaps, denoted by \code{NA}.
#' @param z2 second time series (possibly) with gaps, denoted by \code{NA}.
#' @param gap1 indexes of missing values from \code{z2}, from \code{1:N}, where \code{N = length(z2)}.
#' @param gap2 indexes of missing values from \code{z1}, from \code{1:N}, where \code{N = length(z1)}.
#' @param maxit maximum number of iterations for convergence in interpolation.
#' @param progress logical: should progress be written to screen as iterations proceed?
#' @param sigClip probabilistic significance for choice of line components, dividing series
#     into ``signal'' and ``noise'' (see algorithm for more). Suggested that this be kept
#     above \code{0.95} at a minimum.
#' @param delT the time step delta-t in seconds.
#' @details This function implements the algorithm developed and explained in Chapter 4 of
#   ``Air Pollution and Health: Time Series Tools and Analysis''. 
#' @return zF the final interpolated series.
#' @export
#' @examples library("tsinterp")
#'           data("flux")
#'           z1 <- flux$SagOrig
#'           z1[which(flux$S == FALSE)] <- NA
#'           z2 <- flux$PentOrig
#'           # Unfortunately, not fast enough to run for CRAN checks
#'           sagInt <- BiVarInt(z1 = z1, z2 = z2, gap1 = which(flux$S == FALSE), 
#'                              gap2 = NULL, maxit = 3, delT = 86400)
#'@details
#'Additional details...
#'           A list of five elements, including an interpolated series:
#'            zF the final interpolated series.
#'            p the number of iterations.
#'            diffC the difference between the final series and the previous iteration (metric for convergence).
#'            zA a list of interim series, showing each stage of the convergence.
#'            converge logical indicating whether convergence occurred.


BiVarInt <- function(z1, z2, gap1, gap2, maxit, progress=FALSE, sigClip=0.999, delT=1) {

  cat("Iteration 0:  N/A  (")
  sfInit(parallel = TRUE, cpus = 2)
  ########################################################################
  #
  #  Series z1 setup
  #
  gapTrueA <- rep(NA, length(z1)) 
  gapTrueA[-gap1] <- TRUE
  blocksA <- findBlocks(gapTrueA)

  # linearly interpolated; starting point
  z0a <- z1
  zIa <- linInt(z0a, blocksA)

  # parameters
  N <- length(z1)

  # estimate Mt0 and Tt0
  MtPa <- estimateMt(x=zIa, N=N, nw=5, k=8, pMax=2)

  TtTmpa <- estimateTt(x=zIa - MtPa, epsilon=1e-6, dT=delT, nw=5, k=8,
                      sigClip=sigClip, progress=progress)
  freqRet <- attr(TtTmpa, "Frequency")
  if(length(freqRet) > 1 | (length(freqRet)==1 & freqRet[1] != 0)) {
    TtPa <- rowSums(TtTmpa) 
  } else {
    TtPa <- TtTmpa
  }

  freqSaveA <- attr(TtTmpa, "Frequency")

  converge <- FALSE
  while(!converge) {
    cat(".")
    MtJa <- estimateMt(x=zIa-TtPa, N=N, nw=5, k=8, pMax=2)
    TtTmpa <- estimateTt(x=zIa-MtJa, epsilon=1e-6, dT=delT, nw=5, k=8, 
                   sigClip=sigClip, progress=progress, freqIn=freqSaveA)
    freqRet <- attr(TtTmpa, "Frequency")
    if(length(freqRet) > 1 | (length(freqRet)==1 & freqRet[1] != 0)) {
      TtJa <- rowSums(TtTmpa) 
    } else {
      TtJa <- TtTmpa
    }

    max1 <- max(abs(MtJa - MtPa))
    max2 <- max(abs(TtJa - TtPa))

    if(max(max1,max2) < 1e-4) {
      converge <- TRUE
      MtFa <- MtJa
      TtFa <- TtJa
    } else {
    #  cat(paste(max(max1,max2), "  ", sep=""))
      MtPa <- MtJa
      TtPa <- TtJa
    }
  } # internal M_t / T_t loop
  cat(") ")

  # have initial M_t and T_t estimates for Series 1
  Mt0a <- MtFa
  Tt0a <- TtFa

  ########################################################################
  #
  #  Series z2 setup
  #
  if(length(gap2) > 0) {
    gapTrueB <- rep(NA, length(z2)) 
    gapTrueB[-gap2] <- TRUE
    blocksB <- findBlocks(gapTrueB)
  } else {
    gapTrueB <- rep(TRUE, length(z2))
    blocksB <- matrix(data=0, nrow=0, ncol=3)
  }
  # linearly interpolated; starting point
  z0b <- z2
  if(length(blocksB[, 1]) > 0) {
    zIb <- linInt(z0b, blocksB)
  } else {
    zIb <- z0b
  }

  # parameters
  N <- length(z2)

  # estimate Mt0 and Tt0
  MtPb <- estimateMt(x=zIb, N=N, nw=5, k=8, pMax=2)

  TtTmpb <- estimateTt(x=zIb - MtPb, epsilon=1e-6, dT=delT, nw=5, k=8,
                      sigClip=sigClip, progress=progress)
  freqRet <- attr(TtTmpb, "Frequency")
  if(length(freqRet) > 1 | (length(freqRet)==1 & freqRet[1] != 0)) {
    TtPb <- rowSums(TtTmpb) 
  } else {
    TtPb <- TtTmpb
  }

  freqSaveB <- attr(TtTmpb, "Frequency")
  cat("(")
  converge <- FALSE
  while(!converge) {
    cat(".")
    MtJb <- estimateMt(x=zIb-TtPb, N=N, nw=5, k=8, pMax=2)
    TtTmpb <- estimateTt(x=zIb-MtJb, epsilon=1e-6, dT=delT, nw=5, k=8, 
                   sigClip=sigClip, progress=progress, freqIn=freqSaveB)
    freqRet <- attr(TtTmpb, "Frequency")
    if(length(freqRet) > 1 | (length(freqRet)==1 & freqRet[1] != 0)) {
      TtJb <- rowSums(TtTmpb) 
    } else {
      TtJb <- TtTmpb
    }

    max1 <- max(abs(MtJb - MtPb))
    max2 <- max(abs(TtJb - TtPb))

    if(max(max1,max2) < 1e-4) {
      converge <- TRUE
      MtFb <- MtJb
      TtFb <- TtJb
    } else {
    #  cat(paste(max(max1,max2), "  ", sep=""))
      MtPb <- MtJb
      TtPb <- TtJb
    }
  } # internal M_t / T_t loop
  cat(") \n")

  # have initial M_t and T_t estimates for Series 1
  Mt0b <- MtFb
  Tt0b <- TtFb

  #######################################################################  
  #
  #  Estimate Wt jointly
  #

  zI2a <- zIa - Mt0a - Tt0a
  zI2b <- zIb - Mt0b - Tt0b

  zI2a[gap1] <- 0.0
  zI2b[gap2] <- 0.0
  y <- zI2a

  # find maximum gap length -- only interpolating the first series
  maxGap <- max(blocksA[,3])
  maxlag <- max(length(y)+2, 100)
  neh <- max(2*maxGap, 50)
  clipMax <- max(abs(max(zI2a)), abs(min(zI2a)))

  Wt0a <- estimateWt(zI2a, zI2b, gapTrueA, gapTrueB, dT=delT, blocksA, 
                      neh, maxlag, clipMax) 

  # xd1 = zI2a; xd2 = zI2b; ok1 = gapTrueA; ok2 = gapTrueB; deltat = delT; blocks = blocksA;
  # neh = neh; clipMax = clipMax;

  ########################################################################
  # 
  #  Initial estimates complete -- Tt and Mt for both series, Wt for the first series
  #
  #  Mt0a, Tt0a, Wt0a; Mt0b, Tt0b
  #
  #  Combination gives 'pilot' x estimate, re-use z0a and z0b for this 
  z0a <- z1
  z0a[gap1] <- Mt0a[gap1] + Tt0a[gap1] + Wt0a[gap1]
  z0b <- z2
  z0b[gap2] <- Mt0b[gap2] + Tt0b[gap2]

  zA <- vector("list", maxit+1)
  zA[[1]] <- cbind(z0a, Mt0a, Tt0a, Wt0a, z0b, Mt0b, Tt0b)

  p <- 1
  converge <- FALSE
  cnv <- TRUE
  diffP <- 1e20

  while(!converge) {
    cat(paste("Iteration ", p, ": ", sep=""))
    z1 <- .BiVarInt2(zIa=z0a, zIb=z0b, gap1=gap1, gap2=gap2, blocks=blocksA, delT=delT, sigClip=sigClip,
                       freqSaveA=freqSaveA, freqSaveB=freqSaveB, progress=progress)
    # zIa=z0a; zIb=z0b; 
    diffC <- max(abs(z1[[1]] - z0a)) 

    if(diffC < 1e-3) {
      cat(paste(formatC(diffC, width=6, digits=6, format='f'), "\n", sep=""))
      converge <- TRUE
      zF <- z1[[1]]
    } else if (abs(diffC - diffP) < 1e-5) {
      cat(paste(formatC(diffC, width=6, digits=6, format='f'), "\n", sep=""))
      converge <- TRUE
      zF <- 0.5 * (z1[[1]] + zA[[p]][[1]])
    } else if(diffC - diffP > 0.2*diffP) {  # this is the case where the diffs are diverging ...
      cat(paste(formatC(diffC, width=6, digits=6, format='f'), "\n", sep=""))
      converge <- TRUE 
      cnv <- FALSE
      zF <- z0a          # in this case, the current interpolation is 'worse' than the 
                         # previous one, we assume ... so return the previous
    } else {
      diffP <- diffC
      cat(paste(formatC(diffC, width=6, digits=6, format='f'), "\n", sep=""))
      p <- p+1
      if(p > maxit) {
        converge <- TRUE
        cnv <- FALSE
      }
      z0a <- z1[[1]]
      z0b <- z1[[5]]
      zF <- z1[[1]]
    }
    zA[[p]] <- z1
  }
  sfStop()
  if(cnv) {
    return(list(zF, p, diffC, zA, converge=TRUE))
  } else {
    return(list(zF, p-1, diffC, zA, converge=FALSE))
  }
}

################################################################################
#
#   .BiVarInt2
#
#   Assumes data has been pre-interpolated at least once
#   and that function is in iterative loop
#
################################################################################
.BiVarInt2 <- function(zIa, zIb, gap1, gap2, blocks, delT, sigClip, freqSaveA, freqSaveB, progress) {

  # setup parameters
  gapTrueA <- rep(NA, length(zIa)) 
  gapTrueA[-gap1] <- TRUE
  N <- length(zIa)

  # estimate Mt0 and Tt0
  MtPa <- estimateMt(x=zIa, N=N, nw=5, k=8, pMax=2)

  TtTmpa <- estimateTt(x=zIa-MtPa, epsilon=1e-6, dT=delT, nw=5, k=8, 
                   sigClip=sigClip, progress=progress, freqIn=freqSaveA)
  freqRet <- attr(TtTmpa, "Frequency")
  if(length(freqRet) > 1 | (length(freqRet)==1 & freqRet[1] != 0)) {
    TtPa <- rowSums(TtTmpa)
  } else {
    TtPa <- TtTmpa
  }

  converge <- FALSE
  while(!converge) {
    MtJa <- estimateMt(x=zIa-TtPa, N=N, nw=5, k=8, pMax=2)
    TtTmpa <- estimateTt(x=zIa-MtJa, epsilon=1e-6, dT=delT, nw=5, k=8, 
                   sigClip=sigClip, progress=progress, freqIn=freqSaveA)
    freqRet <- attr(TtTmpa, "Frequency")
    if(length(freqRet) > 1 | (length(freqRet)==1 & freqRet[1] != 0)) {
      TtJa <- rowSums(TtTmpa) 
    } else {
      TtJa <- TtTmpa
    }

    max1 <- max(abs(MtJa - MtPa))
    max2 <- max(abs(TtJa - TtPa))

    if(max(max1,max2) < 1e-4) {
      converge <- TRUE
      MtFa <- MtJa
      TtFa <- TtJa
    } else {
    #  cat(paste(max(max1,max2), "  ", sep=""))
      MtPa <- MtJa
      TtPa <- TtJa
    }
  } # internal M_t / T_t loop

  # have M_t and T_t estimates
  Mt0a <- MtFa
  Tt0a <- TtFa

  ######################################################################## 
  #
  # second series, MtT and TtT
  #
  gapTrueB <- rep(TRUE, length(zIb))
  blocksB <- matrix(data=0, nrow=0, ncol=3)
 
  N <- length(zIb)

  # estimate Mt0 and Tt0
  MtPb <- estimateMt(x=zIb, N=N, nw=5, k=8, pMax=2)

  TtTmpb <- estimateTt(x=zIb-MtPb, epsilon=1e-6, dT=delT, nw=5, k=8, 
                   sigClip=sigClip, progress=progress, freqIn=freqSaveB)
  freqRet <- attr(TtTmpb, "Frequency")
  if(length(freqRet) > 1 | (length(freqRet)==1 & freqRet[1] != 0)) {
    TtPb <- rowSums(TtTmpb)
  } else {
    TtPb <- TtTmpb
  }

  converge <- FALSE
  while(!converge) {
    MtJb <- estimateMt(x=zIb-TtPb, N=N, nw=5, k=8, pMax=2)
    TtTmpb <- estimateTt(x=zIb-MtJb, epsilon=1e-6, dT=delT, nw=5, k=8, 
                   sigClip=sigClip, progress=progress, freqIn=freqSaveB)
    freqRet <- attr(TtTmpb, "Frequency")
    if(length(freqRet) > 1 | (length(freqRet)==1 & freqRet[1] != 0)) {
      TtJb <- rowSums(TtTmpb) 
    } else {
      TtJb <- TtTmpb
    }

    max1 <- max(abs(MtJb - MtPb))
    max2 <- max(abs(TtJb - TtPb))

    if(max(max1,max2) < 1e-4) {
      converge <- TRUE
      MtFb <- MtJb
      TtFb <- TtJb
    } else {
    #  cat(paste(max(max1,max2), "  ", sep=""))
      MtPb <- MtJb
      TtPb <- TtJb
    }
  } # internal M_t / T_t loop

  # have M_t and T_t estimates
  Mt0b <- MtFb
  Tt0b <- TtFb

  #######################################################################  
  #
  #  Estimate Wt jointly
  #

  zI2a <- zIa - Mt0a - Tt0a
  zI2b <- zIb - Mt0b - Tt0b

  zI2a[gap1] <- 0.0
  zI2b[gap2] <- 0.0
  y <- zI2a

  # find maximum gap length -- only interpolating the first series
  maxGap <- max(blocks[, 3])
  maxlag <- max(length(y)+2, 100)
  neh <- max(2*maxGap, 50)
  clipMax <- max(abs(max(zI2a)), abs(min(zI2a)))

  Wt0a <- estimateWt(zI2a, zI2b, gapTrueA, gapTrueB, dT=delT, blocks, 
                      neh, maxlag, clipMax) 

  zIa[gap1] <- Mt0a[gap1] + Tt0a[gap1] + Wt0a[gap1]
  zIb[gap2] <- Mt0b[gap2] + Tt0b[gap2]
  return(list(zIa, Mt0a, Tt0a, Wt0a, zIb, Mt0b, Tt0b))
}



#######################################################################
#
#    mwXSwiener
#
#    Cross-Spectral Bivariate Interpolator
#
#    Takes pre-interpolated series, makes no assumptions about xd2
#    and interpolates missing values in xd1. Can be called sequentially
#    to gapfill xd2 by simply flipping the two series and recomputing 
#    the ok/ACVF/CCVF series.
#
#######################################################################

mwXSwiener <- function(xd1, xd2, ok1, ok2, R11, R12, R21, R22) {

  # establish tag arrays for both series
  n1 <- length(xd1)
  n2 <- length(xd2)

  if(n1 != n2) { stop("Series 1 and Series 2 are of different lengths.") }

  tag1 <- which(!is.na(ok1))
  neq1 <- length(tag1)

  tag2 <- which(!is.na(ok2))
  neq2 <- length(tag2)

  #print(tag1)
  #print(tag2)
  # Block Matrix A
  Amat <- matrix(data=0, nrow=neq1, ncol=neq1)
  for(p in 1:neq1) {
    idx <- tag1[p] - tag1
    Amat[p, ] <- R11[idx + attr(R11, "ZeroF")]
  }

  # Block Matrix B
  Bmat <- matrix(data=0, nrow=neq1, ncol=neq2)
  for(p in 1:neq1) {
    idx <- tag1[p] - tag2
    Bmat[p, ] <- R12[idx + attr(R12, "ZeroF")]
  }

  # Block Matrix C
  Cmat <- matrix(data=0, nrow=neq2, ncol=neq1)
  for(p in 1:neq2) {
    idx <- tag2[p] - tag1
    Cmat[p, ] <- R21[idx + attr(R21, "ZeroF")] 
  }

  # Block Matrix D
  Dmat <- matrix(data=0, nrow=neq2, ncol=neq2) 
  for(p in 1:neq2) {
    idx <- tag2[p] - tag2
    Dmat[p, ] <- R22[idx + attr(R22, "ZeroF")]
  }

  Zmat <- rbind(cbind(Amat, Bmat), cbind(Cmat, Dmat))

  yd1 <- xd1
  yd2 <- xd2
  # Can be done in parallel ?
  for(m in 1:n1) {
    if(is.na(ok1[m])) {
      idx1 <- m - tag1
      wts1 <- R11[idx1 + attr(R11, "ZeroF")]
      idx2 <- m - tag2
      wts2 <- R22[idx2 + attr(R22, "ZeroF")]
      bMat <- c(wts1, wts2)
      wts <- qr.solve(Zmat, bMat) 
      # colnames(wts) <- c(idx1, idx2)
      yd1[m] <- wts[1:neq1] %*% yd1[tag1] + wts[(neq1+1):(neq1+neq2)] %*% yd2[tag2]
    }
  }
  yd1
}



