% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yield.cleaning.R
\name{clean.overlaps}
\alias{clean.overlaps}
\title{Correct for harvester overlaps}
\usage{
clean.overlaps(yield, swath.width, plot = FALSE)
}
\arguments{
\item{yield}{Input yield SpatialPointsDataFrame with column 'yield'}

\item{swath.width}{Swath widht of harvester}

\item{plot}{TRUE/FALSE Plot yield data and detected errors}

\item{col}{colors for plotting}
}
\value{
Returns a SpatialPointsDataFrame with added column 'yield.adj' containing error 
corrected yield.
}
\description{
Correct harvester overlaps by creating polygons for each harvest regions, adjusting them
to remove for areas already harvested and scaling yield accprdingly.
}
\author{
Fiona Evans
}
