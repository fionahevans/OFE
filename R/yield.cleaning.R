###
### Functions for cleaning raw yield data
### 
###
### Expects columns in input SpatialPointsDataFrame: yield, speed, width, distance
###


#' Correct for harvester overlaps
#' 
#' Correct harvester overlaps by creating polygons for each harvest regions, adjusting them
#' to remove for areas already harvested and scaling yield accprdingly.
#'
#' @param yield Input yield SpatialPointsDataFrame with column 'yield' 
#' @param swath.width Swath widht of harvester
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column 'yield.adj' containing error 
#' corrected yield.
#' @export
clean.overlaps <- function(
  yield,
  swath.width,
  plot = FALSE)
{
  if (plot) plot(yield, pch='.', col="gray")
  
  yield$yield.adj <- rep(NA, nrow(yield))
  
  for (i in 2:nrow(yield)) {  
    
    p1 <- coordinates(yield)[i-1, ]
    p2 <- coordinates(yield[i, ])
    
    dist <- sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    theta <- atan2(p2[2] - p1[2], p2[1] - p1[1]) 
    
    # Create rectangle at (0,0) with angle pi/2
    x1 <- c(-dist/2, swath.width/2)
    x2 <- c(dist/2, swath.width/2)
    x3 <- c(dist/2, -swath.width/2)
    x4 <- c(-dist/2, -swath.width/2)
    
    # Rotate 
    x1 <- my.rotate(x1, theta, p2)
    x2 <- my.rotate(x2, theta, p2)
    x3 <- my.rotate(x3, theta, p2)
    x4 <- my.rotate(x4, theta, p2)
    
    if (plot) {
      lines(rbind(x1, x2, x3, x4, x1))
      polygon(rbind(x1, x2, x3, x4, x1), 
              col=color.of(yield$yield[i], range=quantile(yield$yield, c(0, 0.99))))
    }
    
    # Make polygon 
    ply <- mapview::coords2Polygons(rbind(x1, x2, x3, x4, x1), ID = i, 
                                    data = yield@data[i, c("yield", "yield.adj")])
    # ply@polygons[[1]]@Polygons[[1]]@coords
    crs(ply) <- crs(yield)
    
    if (i == 2) polys <- ply
    
    if (i > 2 && rgeos::gIsValid(ply)) {
      
      # Get intersection of current ply with polys
      z <- intersect(polys, ply)
      
      if (!is.null(z)) {
        
        # Remove intersection 
        p <- rgeos::gDifference(ply, z)
        
        # If p is invalid, try buffering
        if (!rgeos::gIsValid(p)) p <- rgeos::gBuffer(p, byid=TRUE, width=0)
        
        if (!is.null(p) && rgeos::gIsValid(p)) {
          
          # Convert back to SpatialPolygonsDataFrame
          row.names(p) <- row.names(ply)
          ply1 <- SpatialPolygonsDataFrame(p, data = ply@data)
          if (plot) plot(ply1, add=T, lwd=2)
          
          # Add to polys
          polys <- rbind(polys, ply1)
          
          # Rescale yield
          yield$yield.adj[i] <- yield$yield[i] * area(ply1) / area(ply)
        }
        
      }
      if (is.null(z)) {
        # Add to polys
        polys <- rbind(polys, ply)
      }
    }
    
  }
  
  yield
}


#' Correct flow delay errors
#' 
#' Correct errors due to grain flow delay. Method adapted from
#' Lee, D. H., et al. (2012). "Automated yield map delay identification using phase correlation methodology." 
#' Transactions of the ASABE 55(3): 743-752.
#'
#' @param yield Input yield SpatialPointsDataFrame with column 'yield'
#' @param delays Vector of delays to test
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column containing error flags.
#' @export
clean.flow <- function(
  yield, 
  delays = c(0:10), 
  plot = FALSE, 
  col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")) 
{
  
  my.lag <- function (x, lag) {
    if (lag == 0) ret <- x
    if (lag > 0) ret <- c(rep(NA, lag), x[1:(length(x)-lag)])
    ret
  }
    
  M <- matrix(nrow=length(delays), ncol=3)
  for (i in 1:length(delays)) {
    print(i)
    # Lag yield
    yield$lag <- my.lag(yield$yield, delays[i])
    # Extimate semivariogram
    v <- variogram(lag ~ 1, data=yield[!is.na(yield$lag), ])
    M[i,] <- as.numeric(v[1, 1:3])
  }
  # Find lag with minimal variance 
  i <- which(M[,3] == min(M[,3]))
  yield$yieldf <- my.lag(yield$yield, delays[i])
  
  if (plot) {
    par(mfrow=c(1,2))
    plot(yield, pch='.', 
         col=color.of(yield$yield, col=col, range=quantile(yield$yieldf, c(0.01, 0.99), na.rm=T)))
    title("Uncorrected")
    plot(yield, pch='.', 
         col=color.of(yield$yieldf, col=col, range=quantile(yield$yieldf, c(0.01, 0.99), na.rm=T)))
    title("Corrected")
  }
  
  cat("Flow delay: ", delays[i], "\n")
  
  yield$yield <- yield$yieldf
  yield[!is.na(yield$yield), ]
}

#' Remove geometry duplicates
#' 
#' Remove geometry duplicates
#'
#' @param yield Input yield SpatialPointsDataFrame with column 'yield'
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column containing error flags.
#' @export
clean.geometry <- function(
  yield, 
  plot = FALSE, 
  col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")) 
{
  indx <- duplicated(coordinates(yield))
  if (plot) {
    plot(yield, pch='.', col=color.of(yield$yield, col=col, range=quantile(yield$yield, c(0.01, 0.99))))
    points(coordinates(yield)[indx, 1], coordinates(yield)[indx, 2])
    title("Geometry duplicates")
  }
  yield$error.geometry <- rep(F, nrow(yield))
  yield$error.geometry[indx] <- T
  cat("Geometry errors: ", sum(yield$error.geometry), "\n")

  yield
}

#' Remove points near paddock boundary
#' 
#' Remove points near paddock boundary
#'
#' @param yield Input yield SpatialPointsDataFrame with column 'yield'
#' @param distance Distance threshold (units as per CRS of yield) from paddock boundary 
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column containing error flags.
#' @export
clean.boundary <- function(
  yield, 
  distance = 10,
  plot = FALSE, 
  col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")) 
{
  hull <- SpatialPoints(concaveman(coordinates(yield), concavity=1.5))
  proj4string(hull) <- proj4string(yield)
  
  hullpts <- as.data.frame(hull)
  hullpts <- hullpts[c(1:nrow(hullpts)) %% 5 == 0, ] # Remove every fifth point 
  cat("Processing boundary... ")
  
  d <- distances(coordinates(yield)[, 1], coordinates(yield)[, 2], 
                 hullpts[, 1], hullpts[, 2])
  indx <- apply(d < distance, 1, sum) 
  indx <- indx > 0
 
  if (plot) {
    plot(yield, pch='.', col=color.of(yield$yield, col=col, range=quantile(yield$yield, c(0.01, 0.99))))
    points(coordinates(yield)[indx, 1], coordinates(yield)[indx, 2])
    title("Boundary errors")
  }
  yield$error.boundary <- rep(F, nrow(yield))
  yield$error.boundary[indx] <- T
  cat("Boundary errors: ", sum(yield$error.boundary), "\n")
  
  yield
}

#' Remove velocity outliers and rapid changes in velocity
#' 
#' Remove velocity outliers and rapid changes in velocity
#'
#' @param yield Input yield SpatialPointsDataFrame with columns 'yield' and 'speed'
#' @param quantiles Minimim and maximum velocity quantiles 
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column containing error flags.
#' @export
clean.velocity <- function(
  yield, 
  quantiles = c(0.01, 0.99), 
  plot = FALSE, 
  col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")) 
{
  indx <- yield$speed <= 
    quantile(yield$speed, quantiles[1]) | yield$speed >= quantile(yield$speed, quantiles[2])
  vdiff <- diff(yield$speed)
  indx1 <- c(FALSE, vdiff <= quantile(vdiff, quantiles[1]) | vdiff >= quantile(vdiff, quantiles[2]))
  indx <- indx | indx1
  if (plot) {
    plot(yield, pch='.', col=color.of(yield$yield, col=col, range=quantile(yield$yield, c(0.01, 0.99))))
    points(coordinates(yield)[indx, 1], coordinates(yield)[indx, 2])
    title("Velocity errors")
  }
  yield$error.velocity <- rep(F, nrow(yield))
  yield$error.velocity[indx] <- T
  cat("Velocity errors: ", sum(yield$error.velocity), "\n")
  
  yield
}

#' Remove harvest overlaps
#' 
#' Remove harvest overlaps ('small' swath widths)
#'
#' @param yield Input yield SpatialPointsDataFrame with columns 'yield' and 'width', 
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column containing error flags.
#' @export
clean.width <- function(
  yield, 
  plot = FALSE, 
  col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")) 
{
  indx <- yield$width <= max(yield$width) * 0.707
  if (plot) {
    plot(yield, pch='.', col=color.of(yield$yield, col=col, range=quantile(yield$yield, c(0.01, 0.99))))
    points(coordinates(yield)[indx, 1], coordinates(yield)[indx, 2])
    title("Harvest overlaps")
  }
  yield$error.width <- rep(F, nrow(yield))
  yield$error.width[indx] <- T
  cat("Width errors: ", sum(yield$error.width), "\n")
  
  yield
}

#' Remove time lags
#' 
#' Remove time lags (small and large distance between points)
#'
#' @param yield Input yield SpatialPointsDataFrame with columns 'yield' and 'distance'
#' @param quantiles Minimim and maximum distance quantiles
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column containing error flags.
#' @export
clean.lags <- function(
  yield, 
  quantiles = c(0.01, 0.99),
  plot = FALSE, 
  col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")) 
{
  indx <- yield$distance <= quantile(yield$distance, quantiles[1]) | 
    yield$distance >= quantile(yield$distance, quantiles[2])
  if (plot) {
    plot(yield, pch='.', col=color.of(yield$yield, col=col, range=quantile(yield$yield, c(0.01, 0.99))))
    points(coordinates(yield)[indx, 1], coordinates(yield)[indx, 2])
    title("Time lags")
  }
  yield$error.lag <- rep(F, nrow(yield))
  yield$error.lag[indx] <- T
  cat("Time lag errors: ", sum(yield$error.lag), "\n")
  
  yield
}

#' Remove yield outliers
#' 
#' Remove yield outliers (small and large yields)
#'
#' @param yield Input yield SpatialPointsDataFrame with column 'yield'
#' @param quantiles Minimim and maximum yield quantiles
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column containing error flags.
#' @export
clean.outliers <- function(
  yield,
  quantiles = c(0.01, 0.99),
  plot = FALSE, 
  col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")) 
{
  indx <- yield$yield <= quantile(yield$yield, quantiles[1]) | 
    yield$yield >= quantile(yield$yield, quantiles[2])
  if (plot) {
    plot(yield, pch='.', col=color.of(yield$yield, col=col, range=quantile(yield$yield, c(0.01, 0.99))))
    points(coordinates(yield)[indx, 1], coordinates(yield)[indx, 2])
    title("Yield outliers")
  }
  yield$error.outlier <- rep(F, nrow(yield))
  yield$error.outlier[indx] <- T
  cat("Yield outlier errors: ", sum(yield$error.outlier), "\n")
  
  yield
}


#' Clean yield data
#' 
#' Remove yield outliers (small and large yields)
#'
#' @param yield Input yield SpatialPointsDataFrame with columns 'yield', 'speed', 'width' and 
#' 'distance'
#' @param boundary.distance Distance from boundary to remove
#' @param quantiles Minimim and maximum yield quantiles
#' @param plot TRUE/FALSE Plot yield data and detected errors
#' @param col colors for plotting
#'
#' @author Fiona Evans
#' 
#' @return Returns a SpatialPointsDataFrame with added column containing error flags.
#' @export
clean <- function(
  yield, 
  boundary.distance = NULL,
  quantiles = c(0.01, 0.99),
  plot = FALSE, 
  col = two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")) 
{
  if (plot) par(mfrow = c(2, 3))
  yield <- OFE::clean.geometry(yield, plot, col)
  if (!is.null(boundary.distance)) yield <- OFE::clean.boundary(yield, distance=boundary.distance, plot, col)
  else yield$error.boundary <- rep(FALSE, nrow(yield))
  yield <- OFE::clean.velocity(yield, quantiles, plot, col)
  yield <- OFE::clean.width(yield, plot, col)
  yield <- OFE::clean.lags(yield, quantiles, plot, col)
  yield <- OFE::clean.outliers(yield, quantiles, plot, col)
  
  indx <- yield$error.geometry | yield$error.boundary | yield$error.velocity |
          yield$error.width | yield$error.lag | yield$error.outlier
  
  cat("Total errors: ", sum(indx), "\n")
  yield <- yield[!indx, ]
  par(mfrow = c(1, 1))
  
  yield
}
