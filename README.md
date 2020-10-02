OFE:  R package containing functions for analysis of data from large on-farm experiments.
====================================================

Functions for cleaning yield monitor data:
* clean: Clean yield data.
* clean.flow: Correct flow delay errors.
* clean.geometry: Remove geometry duplicates.
* clean.lags: Remove timelags.
* clean.outliers: Remove yiled outliers.
* clean.overlaps: Correct for harvester overlaps (EXPERIMENTAL)
* clean.velocity: Remove velocity outliers and rapid changes in velocity
* clean.width: Remove harvest overlaps identifeid by small swath width measurements

Functions for performing ordinary co-kriging for OFE using gstat:
* global.cokrig: Global standardised ordinary co-kriging for OFE using gstat.
* local.cokrig: Local standardised ordinary co-kriging for OFE using gstat.
* cokrig.regression: Performs local linear regressions on co-kriged OFE data.

Functions for plotting data
* color.of: Makes colors for plotting.
* legend.col: Creates a color legend.
* two.colors.old: Older version of 'two.colors' from package 'fields' (less of the middle colour present)
* plot.sim: Plot simulated OFE SpatialPointsDataFrame (as raster).
* plot.cokrig: Plot co-kriged predictions for simulated data.
* plot.cokrig.regression: Plot co-kriged regression results.

Functions for making interpolated raster maps from point data:
* spatial.interp.inla: Int.erpolate to a grid using R-INLA


Data:
* barley.yield: Sample yield monitor data.
* data.linear: Simulated data with global linear response.
* data.linear.zones: Simulated data with linear response that varies for three different zones.
* data.linear.local: 	Simulated data with locally-varying linear response

