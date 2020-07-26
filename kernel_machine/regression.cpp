// This code was written by C.Meixide-Garcia

#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


double KernelG(vec x1,vec x2,float tau) {
  
  int m=x1.n_elem;
  float taum2=1./pow(tau,2);
  vec d;
  
  d=x1-x2;
  
  return exp(-sqrt(dot(d,d))*taum2);
  
}



vec coef(mat X,vec Y,mat W,double lambda,vec mu) {
  
  int i=0;
  int n=X.n_rows;
  int p=X.n_cols;
  vec alphatilde=zeros(n);
  vec beta=zeros(p);
  
  
  mat K=zeros(n,n);
  

  
  alphatilde=inv(K+lambda*eye(n,n))*(W*Y+(eye(n,n)-W)*mu);
  return(alphatilde);}

// [[Rcpp::export]]

double kernel_machine(mat X,vec Xpred,vec coef,vec grid) {
  int n=X.n_rows;
  int i;
  double f=0;
  double Ki=0;
  mat K=zeros(n,n);
  
  
  for(i=0;i<n;i++){
    Ki=KernelG(X.row(i),Xpred,2); //OLLO NON VALE PONHER tau=2, como calculo tau?
    f=f+Ki*coef(i);
  }
  
  
  return(f);
  
}