# https://github.com/CSIRO-Precision-Agriculture/PAT_QGIS_Plugin/blob/master/pat/R-Scripts/_cokrige_Whole_of_Block_Analysis.R
# https://publications.csiro.au/rpr/download?pid=csiro:EP194679&dsid=DS5
# https://usermanual.wiki/Document/PATUserManual.1711714387/html


#' Global standardised ordinary co-kriging for OFE using gstat.
#' 
#' Global standardised ordinary co-kriging for OFE using gstat.
#'
#' @param data SpatialPointsDataFrame containing input data.
#' @param predpts  SpatialPointsDataFrame coordinates to predict at.
#' @param yvar name of varaible to predict.
#' @param xvar name of treament variable.
#' @param model model type, e.g. "Exp", "Sph", "Gau", "Mat" (see gstat help).
#' @param cutoff maxium distance for variogram estimation (see variogram help).
#' @param width width of treatments for calculatign cutoff.
#'
#' @author Fiona Evans
#' 
#' @return Returns SpatialPointsDataFrame containing predictions, variances and covariances. 
#' @export
global.cokrig <- function(data, pred.points, yvar = "yield", xvar = "treatment",
                          model = "Sph", cutoff = NULL, width = NULL){
  
  if (is.null(cutoff) && is.null(width)) {
    stop("Either cutoff or width must be specified.")
  }
  
  # Treatments
  treatments <- sort(unique(data@data[,xvar]))
  nt <- length(treatments)
  nms <- rep(NA, nt)
  
  # Make columns for each treatment
  for (i in 1:nt) {
    nm <- paste0(xvar, treatments[i])
    nms[i] <- nm
    data@data[, nm] <- rep(NA, nrow(data))
    j <- data@data[,xvar] == treatments[i]
    data@data[j, nm] <- data@data[j, yvar]
  }
  
  # Make gstat object
  g <- gstat(NULL, nms[1], form = as.formula(paste(nms[1], '~' ,1)), 
             data = data[data@data[, xvar] == treatments[1],], 
             set = list(nocheck = 1))
  for (i in 2:nt) {
    g <- gstat(g, nms[i], form = as.formula(paste(nms[i], '~' ,1)), 
               data = data[data@data[, xvar] == treatments[i],], 
               set = list(nocheck = 1))
  }
  
  if ( is.null(cutoff) ) {
    cutoff <- width * (nt - 1) + 1
  }
  
  # Estimate variogram
  nBins    <- 20;
  v0 <- variogram(g, cutoff = cutoff, width = cutoff/nBins)
  v <- fit.lmc(v = v0, g = g, vgm(model = model))
  #v$set=list(nocheck=1)
  
  # Predict (hiding gstat print statements)
  invisible(capture.output(p <- predict(v, pred.points)))
  p

}

#' Local standardised ordinary co-kriging for OFE using gstat.
#' 
#' Local standardised ordinary co-kriging for OFE using gstat.
#'
#' @param data SpatialPointsDataFrame containing input data.
#' @param predpts  SpatialPointsDataFrame coordinates to predict at.
#' @param yvar name of varaible to predict.
#' @param xvar name of treament variable
#' @param model model type, e.g. "Exp", "Sph", "Gau", "Mat" (see gstat help).
#' @param max.dist Distance threshold specifying local neighbourhood for around a prediction point,
#' estimating local variogram and making prediction. Must be large enough to inlcude all treatments 
#' at all prediction points.
#'
#' @author Fiona Evans
#' 
#' @return Returns SpatialPointsDataFrame containing predictions, variances and covariances. 
#' @export
local.cokrig <- function(data, pred.points, yvar = "yield", xvar = "treatment",
                          model = "Sph", max.dist = 40){
  
  # Treatments
  treatments <- sort(unique(data@data[,xvar]))
  nt <- length(treatments)
  
  # Distances between data and pred.points
  d.gr <- fields::rdist(coordinates(data), coordinates(pred.points))
  
  # Return object
  ret <- matrix(NA, nrow = ncol(d.gr), 
                ncol = 2 + nt*2 + factorial(nt) / factorial(nt-2)/factorial(2))
  
  # Loop through prediction points, estimate local variogram and predict
  for (i in 1:ncol(d.gr)) {
    #cat(i, "\n")
    
    tmp <- which(d.gr[, i] <= max.dist)
    dat <-  data[tmp, ]
    
    # Check all treatments included 
    n <- table(dat@data[, xvar])
    try( if (any(n < 2) || length(n) < nt) stop("Error: 'max.dist too small."))
    
    z <- as.data.frame(global.cokrig(dat, pred.points[i,], yvar, xvar, model, cutoff = max.dist))[1, ]
    
    ret[i, ] <- as.numeric(z)
  }
  ret <- as.data.frame(ret)
  names(ret) <- names(z)
  coordinates(ret) <- ~x + y
  
  ret
}

#' Plot co-kriged predictions for simulated data.
#' 
#' Plot co-kriged predictions of simulated data.
#'
#' @param data SpatialPointsDataFrame prediction on grid from global.cokrig / local.cokrig 
#' @param rates Vector of treatment names / rates.
#'
#' @author Fiona Evans
#' 
#' 
#' @export
plot.cokrig <- function(data, rates){
  d <- as.data.frame(data)
  gridded(d) <- ~x + y
  d <- stack(d)
  
  nms <- names(d)[grepl("\\pred$", names(d))]
  plot(d[[nms]], nc=4, axes=F, box=F,
       legend.mar=15, 
       axis.args=list(cex.axis=1),
       legend.width=1.2, legend.shrink=0.75,
       col = rev(terrain.colors(99)),
       zlim = c(min(values(d[[nms]])), max(values(d[[nms]]))))
}

#' Performs local linear regressions on co-kriged OFE data.
#' 
#' Performs local linear regressions on co-kriged OFE data. 
#' Fixes the Intercept to the control / zero rate.
#'
#' @param data SpatialPointsDataFrame output (on a grid) from global.cokrig or local.cokrig.
#' @param rates Vector of treatment names (first element should be the control / zero rate).
#'
#' @author Fiona Evans
#' 
#' @return Returns raster stack containing Intercept and Slope. 
#' @export
cokrig.regression <- function(data, rates) {
  data <- as.data.frame(data)
  gridded(data) <- ~x + y
  data <- stack(data)
  
  nms <- names(data)[grepl("\\pred$", names(data))]
  
  # Intercept is zero rate, so just fit slope
  fun <- function(x) { 
    lm(x ~ rates)$coef[2] 
  }
  rate <- raster::calc(data[[nms]] - data[[nms[1]]], fun)
  ret <- stack(data[[nms[1]]], rate)
  names(ret) <- c("Intercept", "Slope")
  ret
}

#' Plot co-kriged regression results.
#' 
#' Plot co-kriged regression results.
#'
#' @param data Raster stack output from cokrig.regression().
#'
#' @author Fiona Evans
#' 
#' @export
plot.cokrig.regression <- function(data){
  plot(data, nc=4, axes=F, box=F,
       legend.mar=15, 
       axis.args=list(cex.axis=1),
       legend.width=1.2, legend.shrink=0.75)
}
