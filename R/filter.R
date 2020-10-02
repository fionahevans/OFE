
#' Filter data to remove low and high outliers
#' 
#' Filter data to removel ow and high outliers
#'
#' @param x Input data.
#' @param low Value of low quantile (default = 0.01).
#' @param high Value of high quantile (default = 0.99).
#'
#' @author Fiona Evans
#' 
#' @return Returns a vector of station ids, in order of increasing distance.
#' @export
qfilter <- function(x, low=0.01, high=0.99) {

  q1 <- quantile(x, low)
  q2 <- quantile(x, high)
  
  x[x <= q1 | x >= q2] <- NA
  x
}
