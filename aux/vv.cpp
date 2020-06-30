#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>



// This code was written by C.Meixide-Garcia

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat vv(arma::mat w,arma::uvec indexw) {
  arma::mat wlittle=w.submat(indexw,indexw);
  return wlittle;
}
  
  
  