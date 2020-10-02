agric:  R package containing functions for accessing & handling agricultural data.
====================================================

Functions:
* SMSyieldCSV: Read yield data from CSV file exported from SMS Advanced software
* swQuery: Query the SkyWatch API for satellite imagery and climate/atmospheric datasets
* swDate: Convert time from SkyWatch API query to a Date object
* swDownload: Download Landsat data using the SkyWatch API and save bands as a raster stack
* getLandsat8LL:	Download Landsat data using the SkyWatch API and save bands as a raster stack
* getLandsat8UTM:	Download Landsat data using the SkyWatch API and save bands as a raster stack
* qfilter:
* grid:	Interpolate data to a grid by kriging and mask using the convex hull of the data
* gwrgrid:	Perform geographically weighted regression on gridded data
* gaussianKernel: Return a 2d array of gaussian weights for use in gwrgrid

APIs for accessing SkyWatch functions are slight modificationson functions from R package 'SkyWatchr'. See the SkyWatchr::querySW help for details on usage.

Planned future functionality:
* extension of gwrgrid to allow for moving window application of any specified R function


## Note:
To compile this package, ensure that the follwoing lines are added to the NAMESPACE file:
```
useDynLib(agric)
importFrom(Rcpp, sourceCpp)
```