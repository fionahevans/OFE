
#' Makes colors for plotting
#' 
#' Makes colors for plotting, stretched to cover the range of x
#'
#' @param x Input data.
#' @param col Color map (default is rainbow).
#' @param range Range of x values to limit color to (not required).
#'
#' @author Fiona Evans
#' 
#' @return Returns a vector of colors.
#' @export
color.of <- function (x, col = rainbow(100), range = NULL, na = NULL) 
{
  if (is.null(range)) 
    range <- range(x, na.rm = T)
  this.col <- col
  if (is.null(na)) this.col <- c("#FFFFFF", col)
  n <- length(this.col)
  indx <- round(my.rescale(x, 2, n, mina = range[1], maxa = range[2]))
  indx <- replace.gt(indx, n, n)
  indx <- replace.lt(indx, 2, 2)
  indx <- na.replace(indx, 1)
  this.col[indx]
}

#' Rescale a vector.
#'
#' Linearly rescale a vector to range between minb and maxb.
#'
#' @param a Input vector.
#' @param minb Minimum of output.
#' @param maxb Maximum of output.
#' @param mina Minimum of input.
#' @param maxa Maximum of input.
#'
#' @keywords manip
#' @export
#' @examples
#' a <- c(2:8)
#' my.rescale(a, 1, 10)
#' my.rescale(a, 1, 10, minb=1, maxb=10)
my.rescale <- function(a, minb, maxb, mina=NULL, maxa=NULL) {
  if (is.null(mina)) mina <- min(a, na.rm=T)
  if (is.null(maxa)) maxa <- max(a, na.rm=T)
  minb + (maxb - minb) * (a - mina)/ (maxa - mina)
}

# Replace values in vector x that are greater than val1 with val2
replace.gt <- function (x, val1, val2) {
  x[x > val1] <- val2
  x
}

# Replace values in vector x that are less than val1 with val2  
replace.lt <- function (x, val1, val2) {
  x[x < val1] <- val2
  x
}

# Replace values in vector x that NA with val 
na.replace <- function(x, val) {
  x[is.na(x)] <- val
  x
}


#' Older version from package 'fields' (less of the middle colour present).
#' 
#' Color interpolation between three colors to output a color vector.
#'
#' @param n Length of output vector.
#' @param start Color.
#' @param end Color.
#' @param middle Color.
#'
#' @keywords color
#' @export
two.colors.old <- function (n = 256, start = "darkgreen", end = "red", middle = "white")
{
  n1 <- n/2
  n2 <- n - n1
  col2 <- col2rgb(end)
  col1 <- col2rgb(start)
  mid.col <- col2rgb(middle)
  e1 <- seq(1, 0, , n1)
  e2 <- seq(0, 1, , n2)
  temp <- rbind(e1 * matrix(col1, nrow = n1, ncol = 3, byrow = TRUE) +
                  (1 - e1) * matrix(mid.col, nrow = n1, ncol = 3, byrow = TRUE),
                e2 * matrix(col2, nrow = n1, ncol = 3, byrow = TRUE) +
                  (1 - e2) * matrix(mid.col, nrow = n1, ncol = 3, byrow = TRUE))
  temp <- temp/256
  rgb(temp[, 1], temp[, 2], temp[, 3])
}

#' Create a color legend.
#' 
#' Legend for a continuous color scale in R, from 
#' https://aurelienmadouasse.wordpress.com/2012/01/13/legend-for-a-continuous-color-scale-in-r/
#'
#' @param lev Data used to produce colors in plot.
#' @param col Color map (default is rainbow).
#'
#' @keywords color
#' @export
legend.col <- function(lev, col=rainbow(100), ystretch=0.5){
  
  opar <- par()
  
  n <- length(col)
  
  # A vector of the form c(x1, x2, y1, y2) 
  # giving the extremes of the user coordinates of the plotting region.
  bx <- par("usr")
  
  len <- (bx[4] - bx[3]) * (1 - ystretch/2)
  
  
  bx[3] <- bx[3] + len
  bx[4] <- bx[4] - len
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3]) 
  box.sy <- (bx[4] - bx[3])  / n
  
  xx <- rep(box.cx, each = 2)
  
  # Draw color boxes
  par(xpd = TRUE)
  for(i in 1:n){
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[n-i+1], border = col[n-i+1])
    
  }
  
  # Add labels
  
  #par(new = TRUE)  
  #par("usr" = bx)

  # plot(0, 0, type = "n",
  #      ylim = c(min(lev, na.rm=T), max(lev, na.rm=T)),
  #      yaxt = "n", ylab = "",
  #      xaxt = "n", xlab = "",
  #      frame.plot = FALSE)
  axis(labels = round(c(min(lev, na.rm=T), max(lev, na.rm=T)), 2), at = c(bx[4], bx[3]), 
       side = 4, las = 2, tick = FALSE, line = .25)
  

}
