#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;

void printMat(NumericMatrix x){
  // print value of matrix
  for (int i=0; i<x.nrow(); ++i) {
    for (int j=0; j<x.ncol(); ++j){
      Rprintf(" %f ", x(i, j));
    }
    Rprintf("\n");
  }
  Rprintf("\n");
}


// [[Rcpp::export]]
SEXP distances(NumericVector x1, NumericVector y1, NumericVector x2, NumericVector y2){
  // Calculate pairwise distances between points
  Rcpp::NumericMatrix ret(x1.size(), x2.size());
  
  for (int i=0; i<x1.size(); i++) {
    for (int j=0; j<x2.size(); j++) {
      ret(i, j) = sqrt(pow(x1[i] - x2[j], 2) + pow(y1[i] - y2[j], 2));
    }
  }

  return ret;
}

NumericVector rotate(double x, double y, double theta, double xcentre, double ycentre){
  Rcpp::NumericVector ret(2);
  ret[1] = - (x + cos(theta) - y * sin(theta));
  ret[2] = - (x + sin(theta) - y * cos(theta));
  
  return ret;
}

// [[Rcpp::export]]
SEXP correctOverlaps(NumericVector x, NumericVector y, NumericVector yield, double width){
  // Adjust yields for area covered using swath width and accounting for overlaps with areas
  // already harvested
//   Rcpp::NumericVector area(x.size);
//   Rcpp::NumericVector adjustedYield(x.size);
//   
//   // Loop through points
//   for (int i=1; i<x.size(); i++) {
//     double dist = sqrt(pow(x[i] - x[i-1], 2) + pow(y[i] - y[i-1], 2));
//     double theta = atan2(y[i] - y[i-1], x[i] - x[i-1]);
//   }
//   
// 
}

/*** R  

# source this file first to test


# Test calling distance from C++
x1 <- y1 <- x2 <- c(1:10)
y2 <- rep(5, 10)
m <- matrix(nrow=10, ncol=10)
for (i in 1:10) {
  for (j in 1:10) {
    m[i, j] <- sqrt((x1[i]-x2[j])^2 + (y1[i]-y2[j])^2)
  }
}
m

# Don't forget to add lines to NAMESPACE before building package:
#   useDynLib(OFE)
#   importFrom(Rcpp, sourceCpp)


distances(as.numeric(x1), as.numeric(y1), as.numeric(x2), as.numeric(y2))

*/
