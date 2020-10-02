#' Interpolate to a grid using R-INLA
#' 
#' Interpolate to a grid using R-INLA
#' Expects columns in input SpatialPointsDataFrame: yield, speed, width, distance
#'
#' @param dat Input SpatialPointsDataFrame with Easting-Northing projection CRS("+proj=utm")
#' @param varname Name of col in dat to interpolate
#' @param width Width of harvester / spreader
#' @param crop TRUE/FALSE Crop concave hull
#' @param plot TRUE/FALSE Plot interpolated data 
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a raster.
#' @export
spatial.interp.inla <- function(dat, 
                                varname, 
                                width, 
                                crop=F, 
                                plot=F, 
                                col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")){
  
  max.edge.length <- c(width, width)*2  # twice header width
  offset <- c(width, width)          # header width
  loc <- as.matrix(coordinates(dat))
  mesh <- inla.mesh.2d(loc = loc, max.edge = max.edge.length, 
                       offset = offset, 
                       cutoff = width)
  
  # Define spatial model (matern)
  spde <- inla.spde2.matern(mesh, alpha = 2) # 0 < alpha < 2 smoothing parameter
  
  # Build stacks for estimation and prediction
  data <- inla.spde.make.A(mesh, loc)
  
  # Define the stack, a different index is needed for each call to the f function
  stk <- inla.stack(data = list(y = dat[, varname]), 
                    A = list(data, 1),
                    effects = list(s = 1:mesh$n, 
                                   data.frame(intercept = rep(1, nrow(dat)))),
                    remove.unused = FALSE, tag = "est")
  
  # Spatial interpolation formula
  formula <- y ~ 0 + intercept + f(s, model=spde) 
  
  
  # Fit the model 
  mod <- inla(formula, data=inla.stack.data(stk),
              control.predictor=list(A = inla.stack.A(stk)),
              family = 'gaussian') 
  
  # Prediction to 10m grid 
  roundDown <- function(x) 10*floor(x/10)
  roundUp <- function(x) 10*ceiling(x/10)
  x <- seq(from=roundDown(min(loc[,1])), to = roundUp(max(loc[,1])), by=10)
  y <- seq(from=roundDown(min(loc[,2])), to = roundUp(max(loc[,2])), by=10)
  grid <- inla.mesh.projector(mesh, xlim=range(x), ylim=range(y), dims=c(length(x), length(y)))
  
  
  pred <- inla.mesh.project(grid, mod$summary.random$s$mean) +
    mod$summary.fixed["intercept", "mean"]
  
  # Rotate and convert to raster
  rot <- t(pred)
  for (i in 1:ncol(rot)) rot[,i] <- rev(rot[,i])
  predr <- raster(rot, xmn=min(x), xmx=max(x), ymn=min(y), ymx=max(y), crs=proj4string(dat))
  
  # Crop areas outside of concave hull
  if (crop) {
    hull <- coords2Polygons(concaveman(coordinates(dat), concavity=1.5), ID="Boundary")
    proj4string(hull) <- proj4string(dat)
    predr <- mask(predr, hull)
  }
  
  if (plot) plot(predr, col=col)
  predr
}