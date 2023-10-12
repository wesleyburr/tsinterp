#' Title estimateMt
#'
#' @param x data?
#' @param N parameter for fitting quadratic
#' @param nw parameter for fitting quadratic
#' @param k parameter for fitting quadratic
#' @param pMax parameter for fitting quadratic
#'
#' @return phat
#' @export
#'
#' @examples
#' 
#' 
estimateMt <- function(x, N, nw, k, pMax) {
  V <- dpss(n=N, nw=5, k=8)$v
  test <- dpssap(V,pMax) # fit quadratic
  U <- test[[1]]
  R <- test[[2]]
  Y <- t(V) %*% x
  a <- t(U) %*% Y
  r <- Y - U %*% a
  xhat <- V %*% r + R %*% a
  phat <- R %*% a
  return(phat)
}