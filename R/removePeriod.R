#' Removes period from a data for time series interpolation purpose
#'
#' @param xd : data
#' @param f0 : freq of periodicity to remove
#' @param nw : parameters of multitaper
#' @param k : parameters of multitaper
#' @param deltaT : parameter of xd
#' @param warn 
#' @param prec : starting precision for finding a good nFFT for removal
#' @param sigClip 
#'
#' @return inv
#' 
#' @export
#' 
#' @importFrom multitaper dpss
#' @importFrom multitaper spec.mtm
#'
#' @examples
#' 
removePeriod <- function(xd, f0, nw, k, deltaT, warn=FALSE, prec=1e-10, sigClip) {
  
  # check to make sure f0 is reasonable, otherwise warn
  N <- length(xd)
  spec.t <- multitaper::spec.mtm(xd,nw=nw,k=k,Ftest=T,plot=F,nFFT=2^(floor(log(N,2))+2),deltat=deltaT)
  idx <- max(which(spec.t$freq < f0))
  if( max(spec.t$mtm$Ftest[idx],spec.t$mtm$Ftest[idx]) < qf(sigClip,2,(2*k-2)) && warn ) {
    warning("Ftest at frequency f0 not significant. Are you sure you want to remove this?")
  }
  
  # early parameter setup, find a nFFT that gives a bin *very* close to f0, or on top of it
  Nyq <- 1/2/deltaT
  nFFT <- -1
  prec.st <- prec
  while( nFFT < 0 ) {
    nFFT <- findPowers(N,f0,Nyq,prec.st)
    prec.st <- prec.st*10
  }
  
  spec <- multitaper::spec.mtm(xd,nw=nw,k=k,returnInternals=T,Ftest=T,plot=F,nFFT=nFFT,maxAdaptiveIterations=0,
                   deltat=deltaT)
  
  # parameter setup
  w <- nw/N/deltaT
  df <- 1/nFFT/deltaT
  neh <- max(10,(floor((2*w)/df+1)))
  f0.idx <- seq(along=spec$freq)[spec$freq == (f0 - min(abs(spec$freq - f0))) | spec$freq == (f0 + min(abs(spec$freq - f0)))]
  
  ##########################################################################
  # 
  #  All spectral window work will require the full spectral array
  # 
  ##########################################################################
  # form spectral windows
  dw <- multitaper::dpss(N,k,5.0)$v*sqrt(deltaT)
  # zero-pad
  dw.z <- rbind(dw,matrix(data=0,nrow=(spec$mtm$nFFT-N),ncol=k))
  # empty window array, nFFT x k
  sw <- matrix(data=0,nrow=spec$mtm$nFFT,ncol=k)
  for(j in 1:k) {
    ft <- fft(dw.z[,j])
    sw[,j] <- c(ft[(spec$mtm$nfreqs+1):spec$mtm$nFFT],ft[1:spec$mtm$nfreqs])
  }
  
  # form estimate of chosen frequency component - takes 0+/- neh from the spectral
  #   window and expands it by multiplying by the CMV at f0
  est <- matrix(data=complex(0,0),nrow=(2*neh+1),ncol=k)
  for(j in 1:k) {
    est[,j] <- spec$mtm$cmv[f0.idx]*(sw[((spec$mtm$nfreqs-1)-neh):((spec$mtm$nfreqs-1)+neh),j])
  }
  
  # subtract from original eigencoefficients
  egn <- spec$mtm$eigenCoefs
  egn <- rbind(Conj(egn[2:spec$mtm$nfreqs,]), egn)
  range <- (f0.idx-neh+spec$mtm$nfreqs) : (f0.idx+neh+spec$mtm$nfreqs)
  if(max(range) > nFFT) {
    # case of folding over the top, i.e. freq too close to Nyq+
    fold <- which(range > nFFT)
    rangeF <- range[fold]
    rangeN <- range[-fold]
    range <- c(1:length(rangeF), rangeN)
  } 
  est2 <- est
  for(j in 1:k) {
    est2[which(range < spec$mtm$nfreqs),j] <- Conj(est2[which(range < spec$mtm$nfreqs),j])
    egn[range,j] <- egn[range,j] - est2[,j]
  }
  
  blank <- matrix(data=0,nrow=nFFT,ncol=1)
  blank[f0.idx] <- spec$mtm$cmv[f0.idx]
  blank[nFFT-f0.idx+2] <- Conj(spec$mtm$cmv[f0.idx])
  inv <- fft(blank,inverse=T)
  inv <- Re(inv)[1:N]
  
  #  cat(paste("Freq: ", spec$freq[f0.idx]," \n",
  #            "Amp : ", sqrt(Mod(spec$mtm$cmv[f0.idx])), "\n",
  #            "Phse: ", Arg(spec$mtm$cmv[f0.idx]), "\n", sep=""))
  return(inv)
}
