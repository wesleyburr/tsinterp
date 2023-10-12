################################################################################
#
#  interpolate
#
#  Sets up gaps, does initial M_t, T_t and W_t estimates, then
#  iterates on interpolate2 until the diff is < 1e-4. 
#
#  ** ADAPT SO PRECISION IS USER-CONTROLLED
#
################################################################################

#' interpolate
#'
#' @param z time series with gaps, denoted by \code{NA}.
#' @param gap indexes of missing values, from \code{1:N}, where \code{N = length(z)}.
#' @param maxit maximum number of iterations for convergence in interpolation.
#' @param progress logical: should progress be written to screen as iterations proceed? 
#' @param sigClip probabilistic significance for choice of line components, dividing series
#'    into ``signal'' and ``noise'' (see algorithm for more). Suggested that this be kept
#'    above \code{0.95} at a minimum.
#' @param delT the time step delta-t in seconds.
#'
#' @return
#' @export
#' Details:   Univariate interpolation of gappy time series. 
#'
#' @examples  
#'    library("tsinterp")
#'    data("flux")
#'    z1 <- flux$SagOrig
#'    z1[which(flux$S == FALSE)] <- NA
#' interpolate(z, gap, maxit = 20, progress=FALSE, sigClip=0.999, delT=1)
#' 
#' 
#' 
#' 
interpolate <- function(z, gap, maxit = 20, progress=FALSE, sigClip=0.999, delT=1) {
  
  stopifnot(is.numeric(delT), delT > 0, 
            is.numeric(sigClip), sigClip > 0, sigClip <= 1.0,
            is.logical(progress),
            is.numeric(maxit), maxit > 0,
            is.numeric(z))
  
  cat("Iteration 0:  N/A  (")
  gapTrue <- rep(NA, length(z)) 
  gapTrue[-gap] <- TRUE
  blocks <- findBlocks(gapTrue)
  
  # linearly interpolated; starting point
  z0 <- z
  zI <- linInt(z0, blocks)
  
  # parameters
  N <- length(z)
  
  # estimate Mt0 and Tt0
  MtP <- estimateMt(x=zI, N=N, nw=5, k=8, pMax=2)
  
  TtTmp <- estimateTt(x=zI - MtP, epsilon=1e-6, dT=delT, nw=5, k=8,
                      sigClip=sigClip, progress=progress)
  freqRet <- attr(TtTmp, "Frequency")
  if(length(freqRet) > 1 | (length(freqRet)==1 && freqRet != 0)) {
    TtP <- rowSums(TtTmp) 
  } else {
    TtP <- TtTmp
  }
  
  freqSave <- attr(TtTmp, "Frequency")
  
  converge <- FALSE
  while(!converge) {
    cat(".")
    MtJ <- estimateMt(x=zI-TtP, N=N, nw=5, k=8, pMax=2)
    TtTmp <- estimateTt(x=zI-MtJ, epsilon=1e-6, dT=delT, nw=5, k=8, 
                        sigClip=sigClip, progress=progress, freqIn=freqSave)
    freqRet <- attr(TtTmp, "Frequency")
    if(length(freqRet) > 1 | (length(freqRet)==1 && freqRet != 0)) {
      TtJ <- rowSums(TtTmp) 
    } else {
      TtJ <- TtTmp
    }
    
    max1 <- max(abs(MtJ - MtP))
    max2 <- max(abs(TtJ - TtP))
    
    if(max(max1,max2) < 1e-4) {
      converge <- TRUE
      MtF <- MtJ
      TtF <- TtJ
    } else {
      #  cat(paste(max(max1,max2), "  ", sep=""))
      MtP <- MtJ
      TtP <- TtJ
    }
  } # internal M_t / T_t loop
  cat(") \n")
  
  # have initial M_t and T_t estimates
  Mt0 <- MtF
  Tt0 <- TtF
  
  # estimate initial W_t
  zI2 <- zI - Mt0 - Tt0
  zI2[gap] <- 0.0
  y <- zI2
  
  # find maximum gap length
  maxGap <- max(blocks[,3])
  neh <- min(5*maxGap, 100)
  clipMax <- max(abs(max(zI2)), abs(min(zI2)))
  # cat(paste("ClipMax = ", clipMax, "\n", sep=""))
  # setup ACV
  spec <- spec.mtm(zI2, nw=5.0, k=8, plot=FALSE, deltat=delT)
  acv <- SpecToACV(spec,maxlag=N)
  # loop on the blocks; NOT AS EFFICIENT?
  # cat(paste("ACV: ", acv[1:4], "\n", sep=""))
  
  for(n in 1:length(blocks[, 1])) { # loop on the blocks
    for(m in blocks[n, 1]:blocks[n, 2]) {
      tag <- ( (blocks[n, 1] - neh):(blocks[n, 2]+neh) )
      tag <- tag[tag %in% (1:N)[-gap]] - m
      bmat <- acv[(abs(tag)+1)]
      neq <- length(tag)
      Amat <- matrix(data=0,nrow=neq,ncol=neq)
      for(q in 1:neq) {
        idx <- abs(tag - tag[q]) 
        Amat[q,] <- acv[idx+1]  # offset because lag-0 == acv[1]
      }
      wts <- qr.solve(Amat, bmat)
      #cat(paste("wts: ", wts, "\n", sep="  "))
      y[m] <- wts %*% y[m+tag]
    }
  }
  Wt0 <- y
  
  ########################################################################
  # 
  #  Initial estimates complete
  #
  #  Mt0, Tt0, Wt0
  #
  #  Combination gives 'pilot' x estimate, re-use z0 for this 
  z0 <- z
  z0[gap] <- Mt0[gap] + Tt0[gap] + Wt0[gap]
  
  zA <- vector("list", maxit+1)
  zA[[1]] <- cbind(z0, Mt0, Tt0, Wt0)
  
  p <- 1
  converge <- FALSE
  cnv <- TRUE
  diffP <- 1e20
  
  while(!converge) {
    cat(paste("Iteration ", p, ": ", sep=""))
    z1 <- .interpolate2(zI=z0, gap=gap, blocks=blocks, delT=delT, sigClip=sigClip,
                        freqSave=freqSave, progress=progress)
    diffC <- max(abs(z1[[1]] - z0)) 
    
    if(diffC < 1e-3) {
      cat(paste(formatC(diffC, width=6, digits=6, format='f'), "\n", sep=""))
      converge <- TRUE
      zF <- z1[[1]]
    } else if (abs(diffC - diffP) < 1e-5) {
      cat(paste(formatC(diffC, width=6, digits=6, format='f'), "\n", sep=""))
      converge <- TRUE
      zF <- 0.5 * (z1[[1]] + zA[[p]][[1]])
    } else if(diffC - diffP > 0.1*diffP) {  # this is the case where the diffs are diverging ...
      cat(paste(formatC(diffC, width=6, digits=6, format='f'), "\n", sep=""))
      converge <- TRUE 
      cnv <- FALSE
      zF <- z0           # in this case, the current interpolation is 'worse' than the 
      # previous one, we assume ... so return the previous
    } else {
      diffP <- diffC
      cat(paste(formatC(diffC, width=6, digits=6, format='f'), "\n", sep=""))
      p <- p+1
      if(p > maxit) {
        converge <- TRUE
        cnv <- FALSE
      }
      z0 <- z1[[1]]
      zF <- z1[[1]]
    }
    zA[[p]] <- z1
  }
  
  if(cnv) {
    return(list(zF, p, diffC, zA, converge=TRUE))
  } else {
    return(list(zF, p-1, diffC, zA, converge=FALSE))
  }
}


################################################################################
#
#   .interpolate2
#
#   Assumes data has been pre-interpolated at least once
#   and that function is in iterative loop
#
################################################################################


#' interpolate2
#'
#' @param zI 
#' @param gap 
#' @param blocks 
#' @param delT 
#' @param sigClip 
#' @param freqSave 
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples
.interpolate2 <- function(zI, gap, blocks, delT, sigClip, freqSave, progress) {
  
  # setup parameters
  gapTrue <- rep(NA, length(zI)) 
  gapTrue[-gap] <- TRUE
  N <- length(zI)
  
  # estimate Mt0 and Tt0
  MtP <- estimateMt(x=zI, N=N, nw=5, k=8, pMax=2)
  
  TtTmp <- estimateTt(x=zI-MtP, epsilon=1e-6, dT=delT, nw=5, k=8, 
                      sigClip=sigClip, progress=progress, freqIn=freqSave)
  freqRet <- attr(TtTmp, "Frequency")
  if(length(freqRet) > 1 | (length(freqRet)==1 && freqRet != 0)) {
    TtP <- rowSums(TtTmp) 
  } else {
    TtP <- TtTmp
  }
  
  converge <- FALSE
  while(!converge) {
    MtJ <- estimateMt(x=zI-TtP, N=N, nw=5, k=8, pMax=2)
    TtTmp <- estimateTt(x=zI-MtJ, epsilon=1e-6, dT=delT, nw=5, k=8, 
                        sigClip=sigClip, progress=progress, freqIn=freqSave)
    freqRet <- attr(TtTmp, "Frequency")
    if(length(freqRet) > 1 | (length(freqRet)==1 && freqRet != 0)) {
      TtJ <- rowSums(TtTmp) 
    } else {
      TtJ <- TtTmp
    }
    
    max1 <- max(abs(MtJ - MtP))
    max2 <- max(abs(TtJ - TtP))
    
    if(max(max1,max2) < 1e-4) {
      converge <- TRUE
      MtF <- MtJ
      TtF <- TtJ
    } else {
      #  cat(paste(max(max1,max2), "  ", sep=""))
      MtP <- MtJ
      TtP <- TtJ
    }
  } # internal M_t / T_t loop
  
  # have M_t and T_t estimates
  Mt0 <- MtF
  Tt0 <- TtF
  
  # estimate W_t
  zI2 <- zI - Mt0 - Tt0
  zI2[gap] <- 0.0
  
  # find maximum gap length
  maxGap <- max(blocks[,3])
  neh <- min(5*maxGap, 100)
  clipMax <- max(abs(max(zI2)), abs(min(zI2)))
  
  # form b matrix
  y <- zI2
  
  # setup ACV
  spec <- spec.mtm(zI2, nw=5.0, k=8, plot=FALSE, deltat=delT)
  acv <- SpecToACV(spec,maxlag=N)
  
  for(n in 1:length(blocks[, 1])) { # loop on the blocks
    for(m in blocks[n, 1]:blocks[n, 2]) {
      tag <- ( (blocks[n, 1] - neh):(blocks[n, 2]+neh) )
      tag <- tag[tag %in% (1:N)[-gap]] - m
      bmat <- acv[(abs(tag)+1)]
      neq <- length(tag)
      Amat <- matrix(data=0,nrow=neq,ncol=neq)
      for(q in 1:neq) {
        idx <- abs(tag - tag[q]) 
        Amat[q,] <- acv[idx+1]  # offset because lag-0 == acv[1]
      }
      wts <- qr.solve(Amat, bmat)
      y[m] <- wts %*% y[m+tag]
    }
  }
  Wt0 <- y
  
  z1 <- zI
  z1[gap] <- Mt0[gap] + Tt0[gap] + Wt0[gap]
  return(list(z1, Mt0, Tt0, Wt0))
}