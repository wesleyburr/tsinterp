#' Compute cross-covariance sequence from two spectral objects.
#'
#' This function computes the cross-covariance sequence from two spectral objects.
# It requires that the spectral objects were computed with the `returnInternals` argument set to TRUE.
# The two spectral objects must have the same parameters for successful computation.
# The function can also limit the maximum lag for the cross-covariance sequence.
# 
#' @param sp1 The first spectral object.
#' @param sp2 The second spectral object.
#' @param maxlag Maximum lag for the cross-covariance sequence (optional).
# 
#' @return A cross-covariance sequence.
# 
#' @details
#' The function computes the cross-covariance sequence by utilizing the spectral properties
#' of two input spectral objects. The `maxlag` parameter controls the maximum lag
#' for the cross-covariance sequence. Ensure that the spectral objects are computed with
#' the `returnInternals` argument set to TRUE and have consistent parameters.
# 
#' @examples
#' sp1 <- spec.mtm(x1, nw = 4, k = 5, returnInternals = TRUE)
#' sp2 <- spec.mtm(x2, nw = 4, k = 5, returnInternals = TRUE)
#' spec.mtm.cross(sp1, sp2)
#'
"spec.mtm.cross" <- function(sp1, sp2, maxlag = NULL) {
  # ... (rest of your function code)
}
