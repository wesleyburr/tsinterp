#' findBlocks

#' Find blocks of missing points from a mask of a time series
#' Scans element-wise to find blocks of missing points
#' 
#' @param mask 
#' Usage: findBlocks(mask)
#' @return a matrix of missing points as blocks, size \code{M * 3}. The start-point, end-point and length of each block are recorded row-wise
#' @export 
#' 
#' 
#' @examples library("tsinterp") 
#'           data("flux")
#'           miss <- flux$S 
#'           miss[miss == FALSE] <- NA
#'           blocks <- findBlocks(miss)
#' 

"findBlocks" <- function(mask) {
  Nlen <- length(mask)
  mask <- which(is.na(mask))
  # case: there are missing points
  if(length(mask) > 0) {
    diffs <- mask[-1] - mask[-length(mask)]
    diffs <- which(diffs > 1)
    
    # case: 1 gap only, possibly no gaps
    if(length(diffs)==0) {
      blocks <- matrix(data=0, nrow=1, ncol=3)
      blocks[1, 1:2] <- c(mask[1], mask[length(mask)])
    } else {
      blocks <- matrix(data=0,nrow=length(mask),ncol=3)
      blocks[1, 1] <- mask[1]
      blocks[1, 2] <- mask[diffs[1]]
      k <- 1
      for(j in 1:length(diffs)) {
        k <- k+1
        blocks[k, 1:2] <- c(mask[diffs[j]+1],mask[diffs[j+1]])
      }
      blocks[k,2] <- max(mask)
      blocks <- blocks[1:k, ]
    }
    blocks[,3] <- blocks[,2] - blocks[,1] + 1
    
    # checks to remove start/end of sequence
    if(blocks[1,1]==1) {
      blocks <- blocks[-1, ]
    }
    if(blocks[length(blocks[,1]),2]==Nlen) {
      blocks <- blocks[-length(blocks[,1]), ] 
    }
  } else {
    blocks <- NULL
  }
  blocks
}
