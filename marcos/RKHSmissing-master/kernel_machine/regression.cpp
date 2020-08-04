// This code was written by C.Meixide-Garcia

#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]

//Compute gaussian kernel 
arma::mat Kernelfull(arma::mat x,float tau) {
  
  int m=x.n_rows;
  float taum2=1./pow(tau,2);
  
  
  arma::mat fil_ig=arma::repmat(x.row(0)*trans(x.row(0)),1,m);
  
  for(int i = 1; i < m; i++){
    fil_ig=join_cols(fil_ig,arma::repmat(x.row(i)*trans(x.row(i)),1,m));
  }
  
  arma::mat col_ig=trans(fil_ig);
  arma::mat cross=x*trans(x);
  
  return arma::exp(-(fil_ig+col_ig-2*cross)*taum2);
  
}

// [[Rcpp::export]]

double KernelG(rowvec x1,vec x2,float tau) {
  
  int m=x1.n_elem;
  float taum2=1./pow(tau,2);
  vec d;
  
  d=x1.t()-x2;
  
  return exp(-(dot(d,d)*taum2));
  
}

// [[Rcpp::export]]

vec coef(mat X,vec Y,mat W,double lambda,vec mu,double sigma) {
  
  int i=0;
  int n=X.n_rows;
  int p=X.n_cols;
  vec alphatilde=zeros(n);
  
  // para que quieres beta 
  vec beta=zeros(p);
  
  
  mat K=Kernelfull(X,sigma); //COMO CALCULO O 2
  mat aux1(n,n); 
  vec aux2(n);
  mat aux3(n,n); 
  aux1=K+lambda*eye(n,n);
  aux3= eye(n,n)-W;
  aux2= W*Y;
  aux2= aux2+aux3*mu;
  aux1= inv(aux1);
  aux2= aux1*aux2;
  aux3= eye(n,n)*lambda+W*K;
  //alphatilde= aux2; 
  alphatilde=inv(K+lambda*eye(n,n))*(W*Y+(eye(n,n)-W)*mu);
  //alphatilde=inv(lambda*eye(n,n)+W*K)*((W*Y));
  return(alphatilde);
  }

// [[Rcpp::export]]

double kernel_machine(mat X,vec Xpred,vec coef,double sigma) {
  int n=X.n_rows;
  int i;
  double f=0;
  double Ki=0;

  
  
  for(i=0;i<n;i++){
    Ki=KernelG(X.row(i),Xpred,sigma); //OLLO NON VALE PONHER tau=2, como calculo tau?
    f=f+Ki*coef(i);
  }
  
  
  return(f);
  
}