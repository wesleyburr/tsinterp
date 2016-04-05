\name{BiVarInt}
\alias{BiVarInt}
\title{Bivariate Interpolation}
\description{
    Bivariate interpolation of gappy time series (possibly two), using a second series
    as a baseline for the interpolation. 
}
\usage{
   BiVarInt(z1, z2, gap1, gap2, maxit, progress=FALSE, sigClip=0.999,
   delT=1, parallelMode = FALSE, ncores = 2)
}
\arguments{
   \item{z1}{first time series (possibly) with gaps, denoted by \code{NA}.}
   \item{z2}{second time series (possibly) with gaps, denoted by \code{NA}.}
   \item{gap1}{indexes of missing values from \code{z1}, from \code{1:N}, where \code{N = length(z1)}. }
   \item{gap2}{indexes of missing values from \code{z2}, from \code{1:N}, where \code{N = length(z2)}. }
   \item{maxit}{maximum number of iterations for convergence in interpolation. }
   \item{progress}{logical: should progress be written to screen as iterations proceed?}
   \item{sigClip}{probabilistic significance for choice of line components, dividing series
         into ``signal'' and ``noise'' (see algorithm for more). Suggested that this be kept
         above \code{0.95} at a minimum.}
   \item{delT}{the time step delta-t in seconds.}
   \item{parallelMode}{logical: if TRUE it will use multicore processing where appliable, if FALSE it will use only sequencial computation}
   \item{ncores}{integer: number of cores to be used if parallel is TRUE }
}
\details{
   This function implements the algorithm developed and explained in Chapter 4 of
   ``Air Pollution and Health: Time Series Tools and Analysis''. 
}
\value{
    A list of five elements, including an interpolated series:
    \item{zF}{the final interpolated series.}
    \item{p}{the number of iterations.}
    \item{diffC}{the difference between the final series and the previous iteration (metric for convergence).}
    \item{zA}{a list of interim series, showing each stage of the convergence.}
    \item{converge}{logical indicating whether convergence occurred.}
}
\examples{
   library("tsinterp")
   data("flux")

   z1 <- flux$SagOrig
   z1[which(flux$S == FALSE)] <- NA
   z2 <- flux$PentOrig

   # Unfortunately, not fast enough to run for CRAN checks
   # sagInt <- BiVarInt(z1 = z1, z2 = z2, gap1 = which(flux$S == FALSE), 
   #                    gap2 = NULL, maxit = 3, delT = 86400)
}

