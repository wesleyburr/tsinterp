\name{flux}
\alias{flux}
\title{flux: Daily Soon Solar Flux}
\description{
    Daily noon solar flux (10.7cm) from Sagamore Hill, Massachusetts
    (USAF monitoring station) and from Penticton, British Columbia (DRAO, NRC). 
    The Penticton data is contiguous, while the Sagamore Hill data contains
    many missing points, indicated by \code{-99} flux. 
}
\format{
    A data frame with six variables: year (\code{Yr}), month (\code{Mnt})
    and day (\code{Day}), two flux series (\code{SagOrig} and \code{PentOrig}),
    and an indicator variable noting which elements of \code{SagOrig} are missing.

    This data set was used as the test case in the development of the algorithms
    in this package.
}
\keyword{datasets}
