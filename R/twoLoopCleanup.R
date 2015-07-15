################################################################################
#
#  twoLoopCleanup
#  ==============
#
#  Takes an input series x, with a provided block length blkL (suggested: length(x)/5 to
#  length(x) / 10), and univariate interpolates the series x with a two-pass system.
#
#  The first pass interpolates chunks as-it-can, eventually getting the entire series
#  interpolated in _some_ form. The second pass re-does the process, interpolating
#  each _gap_ individually, using a large data area around the gap for the best 
#  possible univariate interpolation. 
#
#  Overall, the entire process can take some time, so should only be done on 
#  series with few missing points (e.g., the second "baseline" series for a 
#  bivariate procedure). 
#
#  *** NOTE *** 
#  If the procedure is taking a long time, make blkL smaller, and ensure you are using
#  a series with only a few missing points (e.g., less than 10%). 
#
################################################################################
twoLoopCleanup <- function(x, blkL = 100, parallelMode = TRUE, ncores = 2 ) {
    stopifnot(is.numeric(x), length(x) >= blkL, is.numeric(blkL), blkL != 0)
   # if(blkL < length(x) / 20) { warning("[twoLoopCleanup] blkL is very small: consider increasing it. ") }
  #  if(blkL > length(x) / 5) { warning("[twoLoopCleanup] blkL is very large: consider decreasing it. ") }
  
    # Nothing to interpolate!
    if(length(which(is.na(x))) == 0) { return(x); }
    z0 <- z1 <- x
    ok <- as.numeric(!is.na(x))
    blkL <- floor(blkL)   # cast blkL to integer before using it to compute Nloop, JIC
  
    Nloop <- floor(length(z0) / blkL)

    # first pass: see which things line up on the blocks (not many)
    for(j in 1:Nloop) {
      rng <- ((j-1)*blkL+1):(j*blkL)
      z <- z0[rng]
      gap <- which(ok[rng] == FALSE)   # 15 is hard-coded as an offset from the start/end of 
                                     # the available block
      if(length(gap) > 0) {
        if(length(gap) < blkL/5 & min(gap) > 15 & max(gap) < (blkL-15)) {
          cat("Interpolating: ")
          y <- interpolate(z, gap, maxit=20, progress=FALSE, sigClip=0.99, parallelMode = FALSE, ncores = ncores)
          z1[rng] <- y[[1]]
          ok[rng] <- rep(TRUE, blkL)
        }
      }
    }  #  ** essentially, expect very little to get done in this loop
       #   because there's no guarantee that arbitrarily chosen blocks will line
       #   up with where the actual gaps are!

    z2 <- z1
    finish <- FALSE
    # Step through the gaps, trying to find blocks to work with
    while(!finish) {
      gap <- which(ok == FALSE)  
      diff <- gap[-1] - gap[-length(gap)]
      pos <- which(diff > 25)  # hard-coded as a way to find gaps far enough apart to work with
      if(length(pos) > 0) {
          maxL <- gap[pos[1]+1] - 1
          minL <- max(1, gap[1] - 100)

          # because we're working with a subset ... we need to subset the
          # gaps, because they're absolute
          gap <- gap[gap > minL & gap < maxL]
          gap <- gap - minL + 1
          z <- z2[minL:maxL]
          y <- interpolate(z, gap, maxit=20, progress=FALSE, sigClip=0.99, parallelMode = parallelMode, ncores = ncores)
          z2[minL:maxL] <- y[[1]]
          ok[minL:maxL] <- rep(TRUE, (maxL-minL+1))

          if(length(which(ok==FALSE))==0) {
            finish <- TRUE
          }
      } else {
        finish <- TRUE
      }
    }  # at this point, _most_ of the series should be interpolated; what's left is
       # between floor(length(series) / blkL) and length(series), and anything 
       # that was too tightly clustered to deal with
   
    # stuff at the end of the series
    gap <- which(ok==FALSE)  
    if(length(gap) > 0) {
        minL <- min(gap) - 100; if(minL <= 0) { minL <- 1 }
        maxL <- length(z2)
        gap <- gap[gap > minL & gap < maxL]
        gap <- gap - minL + 1
        z <- z2[minL:maxL]
        y <- interpolate(z, gap, maxit=20, progress=FALSE, sigClip=0.999, parallelMode = parallelMode, ncores = ncores)
        z2[minL:maxL] <- y[[1]]
        ok[minL:maxL] <- rep(TRUE, (maxL-minL+1))
    }

    # 
    #  At this point, the series is fully interpolated, probably badly, or
    #  at best only marginally acceptable. Re-do the interpolation now, using
    #  the first pass as a baseline, instead of the missing points
    #

    #######################################################################
    #
    #  Re-interpolate, one block at a time
    #
    #######################################################################
    # z2 is the current interpolation, and ok is all TRUE

    ok <- !is.na(z0)
    okNA <- ok
    okNA[ok == FALSE] <- NA

    blks <- findBlocks(okNA)   # finds the gaps as start-end-length, so they 
    z3 <- z2                   # can be individually analyzed

    cat("\nRe-interpolating properly: \n\n")
    for(j in 1:length(blks[, 1])) {
      cat(paste("Gap #", j, "\n", sep=""))
      neh <- max(blkL, 3*blks[j, 3])
      rng <- (blks[j, 1] - neh):(blks[j, 2] + neh)
      rng <- rng[rng > 0 & rng <= length(z2)]

      cur <- z2[rng]
      curOK <- rep(TRUE, length(ok))
      curOK[blks[j, 1]:blks[j, 2]] <- FALSE
      gap <- which(curOK[rng]==FALSE)

      y <- interpolate(z=z2[rng], gap, maxit=20, progress=FALSE, sigClip=0.99,  parallelMode = parallelMode, ncores = ncores)
      z3[rng] <- y[[1]]
    }

    z3   # return the interpolated version
}

