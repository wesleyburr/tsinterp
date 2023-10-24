#######################################################################
#
#    SpecToACV
#
#    Converts spectral object to ACVF sequence.
#
#######################################################################


#' SpecToACVdual
#'
#' @param spec 
#' @param maxlag 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
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


#' SpecToACV
#'
#' @param spec 
#' @param maxlag 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
"SpecToACV" <- function(spec,maxlag) {
  s <- spec$spec
  dF <- spec$freq[2] 
  x <- matrix(data=0,nrow=(spec$mtm$nfreqs-1)*2,ncol=1)
  #cat("dim X is ", dim(x), "/n");
  x[1:spec$mtm$nfreqs] = s*dF
  x[(spec$mtm$nfreqs+1):length(x)] <- x[(spec$mtm$nfreqs-1):2]
  x <- as.complex(x)
  x <- Re(fft(x,inverse=TRUE))
  x[1:(maxlag+1)]
}
