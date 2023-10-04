#' Discrete Prolate Spheroidal Sequence Associated Polynomials
#'
#' Generates the associated polynomials for a set of discrete prolate spheroidal
#' sequences (dpss). Used for accurate estimation of polynomial trends.
#'
#' @param V a dpss (Slepian) sequence array, as computed by \code{dpss}.
#' @param maxdeg maximum degree of associated polynomials to estimate and return.
#'
#' @details
#' Computes the Discrete Prolate Spheroidal Sequences Associated 
#' Polynomials for given \code{N} and given (pre-computed) matrix \code{V} 
#' of dimension \code{N * K}, concentration \code{NW}. 
#' Takes parameter \code{maxdeg} as maximum degree.
#'
#' Based on the algorithm developed by Thomson (2001). 
#'
#' @references
#' Thomson, D.J (2001)
#' Spectrum estimation and harmonic analysis. \emph{Proceedings of the IEEE}
#' Volume \bold{70}, number 9, pp. 1055--1096.
#'
#' Thomson, D.J. (2001) Inverse Constrained Projection Filters. 
#' \emph{Proc. SPIE 4478}, Wavelets: Applications in Signal and Image Processing IX, 172 
#' (December 5, 2001); doi:10.1117/12.449708
#'
#' @return
#' Returns a list of three elements, of which the second is the associated polynomials,
#' stored as \code{N * (maxdeg + 1)}. 
#'
#' @examples
#' library("tsinterp")
#' 
#' # compute associated polynomials for given dpss
#' dw <- dpss(n = 100, k = 12, nw = 6)$v
#' dwap <- dpssap(V = dw, maxdeg = 3)
#'
#' @export
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
