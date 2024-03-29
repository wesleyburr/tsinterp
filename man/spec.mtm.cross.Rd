% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spec.mtm.cross.R
\name{spec.mtm.cross}
\alias{spec.mtm.cross}
\title{Compute cross-covariance sequence from two spectral objects.}
\usage{
spec.mtm.cross(sp1, sp2, maxlag = NULL)
}
\arguments{
\item{sp1}{The first spectral object.}

\item{sp2}{The second spectral object.}

\item{maxlag}{Maximum lag for the cross-covariance sequence (optional).}
}
\value{
A cross-covariance sequence.
}
\description{
This function computes the cross-covariance sequence from two spectral objects.
It requires that the spectral objects were computed with the `returnInternals` argument set to TRUE.
The two spectral objects must have the same parameters for successful computation.
The function can also limit the maximum lag for the cross-covariance sequence.
}
\details{
The function computes the cross-covariance sequence by utilizing the spectral properties
of two input spectral objects. The `maxlag` parameter controls the maximum lag
for the cross-covariance sequence. Ensure that the spectral objects are computed with
the `returnInternals` argument set to TRUE and have consistent parameters.
}
\examples{
sp1 <- multitaper::spec.mtm(x1, nw = 4, k = 5, returnInternals = TRUE)
sp2 <- multitaper::spec.mtm(x2, nw = 4, k = 5, returnInternals = TRUE)
spec.mtm.cross(sp1, sp2)

}
