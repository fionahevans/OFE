to.raster <- function(data) {
  data.raster <- as.data.frame(data)
  gridded(data.raster) <- ~x + y
  raster::raster(data.raster)
}

sim.to.stack <- function(data) {
  data.raster <- as.data.frame(data)
  gridded(data.raster) <- ~x + y
  raster::stack(data.raster)
}


#' Convert simulated OFE SpatialPointsDataFrame to raster stack.
#' 
#' Convert simulated OFE SpatialPointsDataFrame to raster stack.
#'
#' @param data SpatialPointsDataFrame with columns 'yield.min', 'rate', 'response', 'yield'.
#'
#' @author Fiona Evans
#' 
#' @return Returns a raster stack
#' @export
sim.to.raster <- function(data) {
  
  stack(to.raster(data[, "yield.min"]),
        to.raster(data[, "yield"]),
        to.raster(data[, "rate"]),
        to.raster(data[, "response"]))
}

#' Plot simulated OFE SpatialPointsDataFrame (as raster).
#' 
#' Plot simulated OFE SpatialPointsDataFrame (as raster).
#'
#' @param data SpatialPointsDataFrame with columns 'yield.min', 'rate', 'response', 'yield'.
#'
#' @author Fiona Evans
#' 
#' @examples
#' data(data.linear.local)
#' 
#' plot.sim(data.linear.local)
#' 
#' @export
plot.sim <- function(data) {
  res <- sim.to.raster(data)
  plot(res, axes=F, box=F, nr=1, nc=5, legend.mar=15, 
       axis.args=list(cex.axis=1),
       legend.width=1.2, legend.shrink=0.75)
}



#' Plot results from basic GWR analysis of simulated data.
#' 
#' Plot results from basic GWR analysis of simulated data.
#'
#' @param data SpatialPointsDataFrame (with columns 'yield.min', 'rate', 'response', 'yield') 
#' used as input to GWmodel::gwr.basic
#' @param model Output from GWmodel::gwr.basic
#'
#' @author Fiona Evans
#' 
#' @export
plot.gwr.basic <- function(dat, model) {
  
  dat$Intercept <- model$SDF$Intercept
  dat$Rate <- model$SDF$rate
  dat$Error <- dat$yield - (model$SDF$Intercept + model$SDF$rate * dat$rate)
  
  dat$Intercept_Error <- dat$yield.min - dat$Intercept 
  dat$Rate_Error <- dat$response - dat$Rate
  
  # Convert to raster stack and plot
  tmp <- cbind(as.data.frame(dat), coordinates(dat))
  gridded(tmp) <- ~x + y
  tmp <- tmp[, c("Intercept", "Rate", "Error", "Intercept_Error", "Rate_Error")]
  
  res <- raster::stack(tmp)
  plot(res, axes=F, box=F, nr=1, nc=5, legend.mar=15, 
       axis.args=list(cex.axis=1),
       legend.width=1.2, legend.shrink=0.75)
}

#' Plot results from INLA spatially-varying covariates analysis of simulated data.
#' 
#' Plot results from INLA spatially-varying covariates analysis of simulated data.
#'
#' @param data SpatialPointsDataFrame (with columns 'yield.min', 'rate', 'response', 'yield') 
#' used as input to fit INLA model
#' @param model Output from INLA::inla
#'
#' @author Fiona Evans
#' 
#' @export
plot.inla <- function(dat, model) {
  # Spatial effect
  dat$Spatial_Effect <- as.numeric(data %*% model$summary.random$s$mean)*1000
  
  # Spatially varying rate effect
  dat$Rate_Effect  <-  as.numeric(data %*% model$summary.random$rate$mean)*1000
  
  # Error and error components
  dat$Error <- dat$yield - 
    (dat$Spatial_Effect + dat$Rate_Effect * dat$rate)
  dat$Spatial_Error <- dat$yield.min - dat$Spatial_Effect 
  dat$Rate_Error <- dat$Rate_Effect - dat$response
  
  # Convert to raster stack and plot
  tmp <- cbind(as.data.frame(dat), coordinates(dat))
  gridded(tmp) <- ~x + y
  tmp <- tmp[, c("Spatial_Effect", "Rate_Effect", "Error", "Spatial_Error", "Rate_Error")]
  res <- raster::stack(tmp)
  plot(res, axes=F, box=F, nr=1, nc=5, legend.mar=15, 
       axis.args=list(cex.axis=1),
       legend.width=1.2, legend.shrink=0.75)
}


#' Simulated data with global linear response.
#'
#' Simulated on-farm experiment data with global linear response to fertiliser.
#' 
#'
#' @docType data
#' @usage data(data.linear)
#' @keywords datasets
#' 
#' @format A SpatialPointsDataFrame  with 2400 rows and 4 variables:
#' \describe{
#'   \item{yield.min}{Minimum yield (no fertiliser applied)}
#'   \item{rate}{Rate of fertiliser applied}
#'   \item{response}{Linear rate of response to fertiliser}
#'   \item{yield}{Resulting yield from experiment}
#' }
#'
#' @references Evans et al. 2020 MDPI Agronomy
#' (\href{https://address here})
#'
#' @examples
#' data(data.linear)
#' head(data.linear)
#' 
#' plot.sim(data.linear)
"data.linear"


#' Simulated data with linear response that varies for three different zones.
#'
#' Simulated on-farm experiment data with linear response to fertiliser 
#' that varies for three different zones.
#' 
#'
#' @docType data
#' @usage data(data.linear.zones)
#' @keywords datasets
#' 
#' @format A SpatialPointsDataFrame  with 2400 rows and 4 variables:
#' \describe{
#'   \item{yield.min}{Minimum yield (no fertiliser applied)}
#'   \item{rate}{Rate of fertiliser applied}
#'   \item{response}{Linear rate of response to fertiliser}
#'   \item{yield}{Resulting yield from experiment}
#' }
#'
#' @references Evans et al. 2020 MDPI Agronomy
#' (\href{https://address here})
#'
#' @examples
#' data(data.linear.zones)
#' head(data.linear.zones)
#' 
#' plot.sim(data.linear.zones)
"data.linear.zones"


#' Simulated data with locally-varying linear response.
#'
#' Simulated on-farm experiment data with locally-varying linear response.
#' 
#'
#' @docType data
#' @usage data(data.linear.local)
#' @keywords datasets
#' 
#' @format A SpatialPointsDataFrame  with 2400 rows and 4 variables:
#' \describe{
#'   \item{yield.min}{Minimum yield (no fertiliser applied)}
#'   \item{rate}{Rate of fertiliser applied}
#'   \item{response}{Linear rate of response to fertiliser}
#'   \item{yield}{Resulting yield from experiment}
#' }
#'
#' @references Evans et al. 2020 MDPI Agronomy
#' (\href{https://address here})
#'
#' @examples
#' data(data.linear.local)
#' head(data.linear.local)
#' 
#' plot.sim(data.linear.local)
"data.linear.local"