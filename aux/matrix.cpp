#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>



// This code was written by C.Meixide-Garcia

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

int grande(int n,arma::mat A,arma::mat B) {
  arma::mat C;
  for(int i=0;i<n;i++){
    C=arma::kron(A,B);
  }
  return 0;
}