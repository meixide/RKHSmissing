// This code was written by C.Meixide-Garcia

#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

double logit(vec Xi,vec beta) {
  return(exp(dot(Xi,beta))/(1+exp(dot(Xi,beta))));
}

vec coef(mat X,vec Y,mat W,double lambda) {
  
  int i=0;
  int n=X.n_rows;
  int p=X.n_cols;
  vec alphatilde=zeros(n);
  vec beta=zeros(p);
  vec mu=zeros(n);
  
  mat K=zeros(n,n);
  
  for (i=0;i<n;i++) {
    mu(i)=(2*logit(X.row(i),beta));
  }
  
  alphatilde=inv(K+lambda*eye(n,n))*(W*Y+(eye(n,n)-W)*mu);
  return(alphatilde);}

// [[Rcpp::export]]

vec kernel_machine(mat X,vec coef,vec grid) {
  int n=X.n_rows;
  mat K=zeros(n,n);
  vec f=K*coef;
  return(f);
  
}