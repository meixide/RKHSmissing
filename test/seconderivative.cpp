#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;

// [[Rcpp::export]]


arma::mat distance_mf(arma::mat x) {
  
  int m=x.n_rows;
  
  
  arma::mat fil_ig=arma::repmat(x.row(0)*trans(x.row(0)),1,m);
  
  for(int i = 1; i < m; i++){
    fil_ig=join_cols(fil_ig,arma::repmat(x.row(i)*trans(x.row(i)),1,m));
  }
  
  arma::mat col_ig=trans(fil_ig);
  arma::mat cross=x*trans(x);
  
  return arma::pow(arma::abs(fil_ig+col_ig-2*cross),0.5);
}

//Compute second derivative of objective function
// [[Rcpp::export]]
arma::mat full_Findweight(int n,int p,arma::mat X,arma::vec Y,float tau,int numberclose,arma::mat wcens) {
  
  /*arma::vec Y,*/
    arma::mat w(n,n);
  arma::mat weight(n,n);
  weight.zeros();
  
  w.zeros();
  weight=w;
  int i=0;
  // float cutoff=0;
  //float dij=0;
  // int j=0;
  float invtau2=1./pow(tau,2);
  
  arma::mat d=arma::zeros(n,n);
  d=arma::pow(distance_mf(X),2.);
  for (i=0;i<n;i++) {
    d(i,i)=arma::datum::inf;
    
  }
  
  
  weight=arma::exp(-d*invtau2);
  
  w=weight;
  w=weight%wcens;
  
  
  
  return(w);
  
}

// [[Rcpp::export]]

arma::umat findneigh(int n,int numberclose,arma::mat d) {
  
  arma::uvec neighi(numberclose);
  arma::uvec d1_index(n);
  arma::umat neigh(n,numberclose);
  arma::rowvec d1(n);
  
  int i=0;
  
  
  for (i=0;i<n;i++) {
    
    d1=d.row(i);
    d1_index=arma::sort_index(d1);
    neighi=d1_index.subvec(1,numberclose);
    neigh.row(i)=arma::trans(neighi);
    
    
  }
  
  
  
  return neigh;
  
}

// [[Rcpp::export]]

arma::mat findwneigh(int n,int numberclose,arma::mat w,arma::umat neigh) {
  
  arma::mat wsmall(n,numberclose);
  arma::rowvec wi(n);
  for (int i=0;i<n;i++) {
    
    wi=w.row(i);
    wsmall.row(i)=arma::trans(wi.elem(neigh.row(i)));
  }
  
  return wsmall;
}

// [[Rcpp::export]]

arma::mat KernelG(arma::mat x,float tau) {
  
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


arma::mat Secondderivative(int n,int p,int numberclose,arma::mat x,arma::vec y,arma::mat wsmall,arma::mat K,arma::umat neigh) {
  arma::mat M=arma::zeros(n*p,n*p);
  arma::mat B=arma::zeros(n,n);
  arma::mat A=arma::zeros(p,p);
  arma::vec yij=arma::zeros(p);
  arma::vec ki=arma::zeros(n);
  
  
  int a,b,c,d=0;
  
  for(int i=0;i<n;i++) {
    ki=K.col(i);
    // B=K.col(i)*arma::trans(K.col(i));
    for(int j=0;j<numberclose;j++) {
      // A=arma::trans()*(x.row(i)-x.row(neigh(i,j)));
      yij=arma::trans(x.row(i)-x.row(neigh(i,j)));
      for(a=0;a<p;a++) {
        for(b=0;b<p;b++) {
          for(c=0;c<n;c++) {
            for(d=0;d<n;d++){
              //M+=2*wsmall[n*j+i] *arma::kron(A,B)*(1./(n*(n-1)));
              
              M(a*n+c,b*n+d)+=2*wsmall[n*j+i]*yij(a)*yij(b)*ki(c)*ki(d)*(1./(n*(n-1)));
            }
          }
        }
      }
    }
  }
  return M;
}
