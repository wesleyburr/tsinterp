\name{linInt}
\alias{linInt}
\alias{linearInterpolation}
\title{Linear Interpolation}
\description{
    Linear interpolation (straight-line interpolation) of a time series.
}
\usage{
    linInt(dat,blocks)
}
\arguments{
   \item{dat}{a time series with length at least \code{max(blocks[, ])}. }
   \item{blocks}{a matrix, typically obtained from \code{\link{findBlocks}{findBlocks}},
         of dimension \code{M * 3}. The first column contains the start-point of the 
         block, the second the end-point of the block, and the third the length of the
         block (the number of missing points in it). }
}
\details{
    Implementation of one of the simplest possible algorithms for interpolation:
    linear interpolation. Each gap is replaced by a straight line interpolation from
    the left to right. 
}
\value{
    The linearly interpolated input \code{dat}. 
}
\examples{
    library("tsinterp")
    data("flux")

    # linearly interpolate the Sagamore Hill flux series
    miss <- flux$S
    miss[miss == FALSE] <- NA
    blocks <- findBlocks(miss)

    fluxInt <- linInt(flux$SagOrig, blocks)
}
