\name{findBlocks}
\alias{findBlocks}
\title{Find blocks of missing points from a mask of a time series.}
\description{
    Given a mask for a time series, scan for missing points (denoted by \code{NA}) and
    block out the gaps so that each can be interpolated in turn.
}
\usage{
    findBlocks(mask)
}
\arguments{
   \item{mask}{an array with \code{NA} and not-\code{NA} elements. The \code{NA} elements
       are considered to be missing, and the not-\code{NA} elements are considered to be
       present. The \code{NA} elements are aggregated and recorded as blocks.}
}
\details{
    Scans element-wise to find blocks of missing points.
}
\value{
    Returns a matrix of missing points as blocks, size \code{M * 3}. The 
    start-point, end-point and length of each block are recorded row-wise.
}
\examples{
   library("tsinterp")

   data("flux")
   miss <- flux$S
   miss[miss == FALSE] <- NA
   blocks <- findBlocks(miss)
}

