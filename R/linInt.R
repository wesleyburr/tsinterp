


#' Title linInt
#'
#' @param dat data
#' @param blocks 
#'
#' @return dat
#' @export
#' 
#' @examples
#' 
"linInt" <- function(dat,blocks) {
  nGap <- length(blocks[,1])
  for(j in 1:nGap) {
    dY <- (dat[blocks[j,2]+1] - dat[blocks[j,1]-1])/(blocks[j,3]+1)
    st <- dat[blocks[j,1]-1]
    lt <- dat[blocks[j,2]+1]
    if(dY != 0) {
      fill <- seq(st,lt,dY)
    } else {
      fill <- rep(st,blocks[j,2]-blocks[j,1]+3) 
    }
    fill <- fill[c(-1,-length(fill))]
    dat[blocks[j,1]:blocks[j,2]] <- fill
  }
  dat
}