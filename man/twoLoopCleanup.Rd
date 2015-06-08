\name{twoLoopCleanup}
\alias{twoLoopCleanup}
\title{Univariate Interpolation: Implementation}
\description{
    Univariate interpolation of gappy time series: implementation wrapping \link{interpolate} in some logic.
}
\usage{
   twoLoopCleanup(x, blkL = 100)
}
\arguments{
   \item{x}{time series with gaps, denoted by \code{NA}.}
   \item{blkL}{length of block to use in first-pass interpolation; suggested as \code{length(x)/10}. }
}
\details{
    This function wraps up calls to \code{interpolate}, another function in this package. Suitable for 
    univariate interpolation of a single series with a small percentage of missing points, say 10 or 20 percent.

    Calls to \code{interpolate} in this routine are constrained to have \code{maxit=20},
    \code{progress=FALSE} and \code{sigClip=0.99}. 
}
\value{
    \item{z}{the final interpolated series.}
}
\examples{
   library("tsinterp")
   data("flux")

   z1 <- flux$SagOrig
   z1[which(flux$S == FALSE)] <- NA

   # Unfortunately, not fast enough to run for CRAN checks
   #sagInt <- twoLoopCleanup(z1, blkL = 100)

}
