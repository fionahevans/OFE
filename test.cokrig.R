library(OFE)

data(data.linear.local)

data <- data.linear.local
pred.points <- data

width <- 3
rates <- c(0, 30, 60, 100)
nt <- length(rates)


data.raster <- sim.to.raster(data)
plot.sim(data)

yvar <- "yield"
xvar <- "rate"

# What would GWR bandwidth be?
bw <- width * (nt - 1) / 2

# Distance to include all treatments
# max.dist is Euclidean distance, for this data distance = equals number of pixels
max.dist <- width * (nt - 1) + 1

model <- "Sph"

# Global standardised ordinary co-kriging ----
system.time(
  gc <- global.cokrig(data, pred.points, yvar = "yield", xvar = "rate",
                            model = "Sph", cutoff = max.dist)
)
# user  system elapsed 
# 68.78    0.17   69.50


plot.cokrig(gc, rates=c(0, 30, 60, 100))

# Local / block standardised ordinary co-kriging ----

system.time(
  lc <- local.cokrig(data, pred.points, yvar = "yield", xvar = "rate",
                model = "Sph", max.dist = max.dist)
)
# user  system elapsed 
# 119.53    0.08  120.48 (2 mins)

plot.cokrig(lc, rates=c(0, 30, 60, 100))


# Lareger max.dist
system.time(
  lc2 <- local.cokrig(data, pred.points, yvar = "yield", xvar = "rate",
                     model = "Sph", max.dist = 20)
)
# user  system elapsed 
# 364.15    4.39  380.97 (6.5 mins)

plot.cokrig(lc2, rates=c(0, 30, 60, 100))


### Local regressions ----

cr <- cokrig.regression(lc, rates = c(0, 30, 60, 100)) 
plot.cokrig.regression(cr)

cr2 <- cokrig.regression(lc2, rates = c(0, 30, 60, 100)) 
plot.cokrig.regression(cr2)


# Look at errors in the intercept and response
par(mfrow=c(1,2))
hist(values(cr[["Intercept"]]) - values(data.raster[["yield.min"]]), xlab="", main="max.dist = 10")
hist(values(cr2[["Intercept"]])  - values(data.raster[["yield.min"]]), xlab="", main="max.dist = 20")

mean(abs(values(cr[["Intercept"]])  - values(data.raster[["yield.min"]])))  
mean(abs(values(cr2[["Intercept"]])  - values(data.raster[["yield.min"]]))) # lower intercept error

par(mfrow=c(1,2))
hist(values(cr[["Slope"]]) - values(data.raster[["response"]]), xlab="", main="max.dist = 10")
hist(values(cr2[["Slope"]])  - values(data.raster[["response"]]), xlab="", main="max.dist = 20")

mean(abs(values(cr[["Slope"]])  - values(data.raster[["response"]])))  # lower response error
mean(abs(values(cr2[["Slope"]])  - values(data.raster[["response"]]))) 
