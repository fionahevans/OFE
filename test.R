# roxygen2::roxygenize('.')
#
# Don't forget to add lines to NAMESPACE before building package:
#   useDynLib(OFE)
#   importFrom(Rcpp, sourceCpp)

library(OFE)
data(barley.yield)

col <- two.colors.old(start = "brown", end = "darkgreen", middle = "yellow")
plot(barley.yield, pch='.', col=color.of(barley.yield$yield, col=col))

# Clean (with boundary cleaning)
y2 <- OFE::clean(barley.yield, boundary.distance = 20, plot = F)

par(mfrow=c(1,2))
plot(barley.yield, pch='.', col=color.of(barley.yield$yield, col=col, 
                                         range=quantile(y2$yield, c(0.01, 0.99))))
plot(y2, pch='.', col=color.of(y2$yield, col=col, 
                               range=quantile(y2$yield, c(0.01, 0.99))))


# Clean (omit boundary cleaning if experiment is not near boundary)
y2 <- OFE::clean(barley.yield, plot = F)

# Plot with same colours
par(mfrow=c(1,2))
plot(barley.yield, pch='.', col=color.of(barley.yield$yield, col=col, 
                                  range=quantile(y2$yield, c(0.01, 0.99))))
plot(y2, pch='.', col=color.of(y2$yield, col=col, 
                                  range=quantile(y2$yield, c(0.01, 0.99))))


# Still testing - this is not fast!
y1 <- OFE::clean.flow(barley.yield, delays = c(0:5), plot = T)

par(mfrow=c(1,2))
plot(yield, pch='.', xlim=c(516000, 516200), ylim=c(6494200, 6494400),
     col=color.of(yield$yield, col=col, range=quantile(yield$yieldf, c(0.01, 0.99), na.rm=T)))
title("Uncorrected")
plot(yield, pch='.', xlim=c(516000, 516200), ylim=c(6494200, 6494400),
     col=color.of(yield$yieldf, col=col, range=quantile(yield$yieldf, c(0.01, 0.99), na.rm=T)))
title("Corrected")




