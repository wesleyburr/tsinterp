#' Estimate the sinusoidal components in a time series.
#'
#' This function estimates sinusoidal components in a time series using the multitaper method.
#'
#' @param x The input time series.
#' @param epsilon A tolerance for peak estimation convergence.
#' @param dT The time interval between data points.
#' @param nw The time-half bandwidth parameter for multitaper.
#' @param k The number of tapers.
#' @param sigClip Significance level for peak detection.
#' @param progress If TRUE, show progress information.
#' @param freqIn A pre-defined set of frequencies (optional).
#'
#' @return A matrix containing sinusoidal components.
#'
#' @details
#' The function follows these algorithmic steps:
#' 1. Spectrum/F-test pilot estimate using the multitaper method.
#' 2. Estimation of significant peaks based on the significance level (sigClip).
#' 3. Refinement of peak locations by eliminating duplicates and optimizing peak frequencies.
#' 4. Estimation of phase and amplitude by inverting the spectrum (line component removal).
#'
#' @examples
#' estimateTt(x, epsilon = 1e-10, dT = 0.01, nw = 4, k = 5, sigClip = 0.05)
#'
estimateTt <- function(x, epsilon, dT, nw, k, sigClip, progress = FALSE, freqIn = NULL) {
  # ... (rest of your function code)
}
