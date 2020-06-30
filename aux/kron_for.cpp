#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

#include <iostream> 
using namespace std;

// This code was written by C.Meixide-Garcia

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat kron_for(arma::mat A,arma::mat B) {
  int p=A.n_cols;
  int n=B.n_cols;
  int a,b,c,d=0;
  Rcpp::Rcout << "p"  << p<< endl; 
  arma::mat M=arma::zeros(n*p,n*p);
    for(a=0;a<p;a++) {
      for(b=0;b<p;b++) {
        for(c=0;c<n;c++) {
          for(d=0;d<n;d++){
            M(a*n+c,b*n+d)=A(a,b)*B(c,d);
            Rcpp::Rcout << "a"  << a<< endl; 
          }
        }
      }
    }
    return M;
}