#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp> //para monitorizar o tempo
#include <progress_bar.hpp>

// [[Rcpp::export]]

vec armakron(vec y,mat A,mat B){
  
  
  // we get a vector of values and a list of matrices from R
  
  
  int nmat =  2;	// how many matrices
  int nall = y.size();
  vec cols(nmat);
  
  // need to map to arma::mat type
  std::vector<mat> list;
  list.push_back( A );
  list.push_back( B );
  // fill in dimensions
  cols(0) = list.at(0).n_cols;
  cols(1) = list.at(1).n_cols;
  
  
  if ( prod(cols) != nall ){
    throw std::logic_error( "armakron: prod(cols(matrices)) not equal length(y)" );
    
  }
  
  // TODO do a sanity check on dimensions of objects!!
  
  // product for first matrix
  vec y0 = y;
  vec y1(nall);
  y1.zeros();
  int n = list.at(0).n_rows;
  int m = nall/n;
  uvec lhs_out = linspace<uvec>(0,n-1,n);
  
  uvec lhs_in = lhs_out * m;
  uvec rhs(n);
  uvec lhs(n);
  
  for (int i=0; i<m; i++){
    lhs = lhs_in + i;
    rhs    = lhs_out + i*n;
    y1.elem( lhs ) = list.at( 0 ) * y0.elem( rhs ) ;
  }
  
  // process all other matrices
  if (nmat > 1){
    for(int imat=1; imat < nmat; imat++){
      y0 = y1;
      n  = list.at(imat).n_rows;
      m  = nall/n;
      lhs_out.resize(n);
      lhs_out = linspace<uvec>(0,n-1,n);
      uvec lhs_in = lhs_out * m;
      rhs.resize(n);
      lhs.resize(n);
      for (int i=0; i<m; i++){
        lhs = lhs_in + i;
        rhs    = lhs_out + i*n;
        y1.elem( lhs ) = list.at( imat ) * y0.elem( rhs ) ;
      }
    }
  }
  return y1;
  
  
}

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
  arma::mat Maux=arma::zeros(n*p,1);
  arma::mat E=arma::zeros(p,p);
  arma::vec canoni=arma::zeros(n*p);
    canoni(0)=1;
for(int k=0;k<n;k++) {
  for(int i=0;i<n;i++) {
    
    B=K.col(i)*arma::trans(K.col(i));
    for(int j=0;j<numberclose;j++) {
    E=E+trans(x.row(i)-x.row(neigh(i,j)))*(x.row(i)-x.row(neigh(i,j)))*wsmall(i,j);
    }
      Maux=Maux+armakron(canoni,B,E);
      //M=M+arma::kron(E,B);
    }
  M.col(k)=Maux;
}
  return M*(2./(n*(n-1)));
  //return Maux*(2./(n*(n-1)));
}

// [[Rcpp::export]]

arma::mat Secondderivativevello(int n,int p,int numberclose,arma::mat x,arma::vec y,arma::mat wsmall,arma::mat K,arma::umat neigh) {
  arma::mat M=arma::zeros(n*p,n*p);
  arma::mat B=arma::zeros(n,n);
  arma::mat A=arma::zeros(p,p);
  
  for(int i=0;i<n;i++) {
    
    B=K.col(i)*arma::trans(K.col(i));
    for(int j=0;j<numberclose;j++) {
      A=arma::trans(x.row(i)-x.row(neigh(i,j)))*(x.row(i)-x.row(neigh(i,j)));
      M=M+2*wsmall(i,j)*arma::kron(A,B)*(1./(n*(n-1)));
    }
  }
  return M;
}
// [[Rcpp::export]]

arma::mat Secondderivativeprobas(int n,int p,int numberclose,arma::mat x,arma::vec y,arma::mat wsmall,arma::mat K,arma::umat neigh,
                                 bool display_progress=true) {
  arma::mat M=arma::zeros(n*p,n*p);
  arma::mat B=arma::zeros(n,n);
  arma::mat A=arma::zeros(p,p);
  int a=0;
  int b=0;
  int c=0;
  int d=0;
  
  Progress q(n, display_progress);
  for(int i=0;i<n;i++) {
    
    B=K.col(i)*arma::trans(K.col(i));
    A=arma::zeros(p,p);
    for(int j=0;j<numberclose;j++) {
      A=A+wsmall(i,j)*arma::trans(x.row(i)-x.row(neigh(i,j)))*(x.row(i)-x.row(neigh(i,j)));
    }
    
    for(a=0;a<p;a++) {
      for(b=a;b<p;b++) {
        for(c=0;c<n;c++) {
          
          if(b==a) {
            for(d=c;d<n;d++){
              
              
              M(a*n+c,b*n+d)+=A(a,b)*B(c,d);
            
          }
          } else {
              for(d=0;d<n;d++){
                
                
                M(a*n+c,b*n+d)+=A(a,b)*B(c,d);
                
              }
              
            }
         
        }
      }
    }
    
      //M=M+arma::kron(A,B);
    q.increment(); 
  }
  return arma::symmatu(M)*(2./(n*(n-1)));
}

// [[Rcpp::export]]
arma::mat Secondderivativearma(int n,int p,int numberclose,arma::mat x,arma::vec y,arma::mat wsmall,arma::mat K,arma::umat neigh,
                               bool display_progress=true) {
  arma::mat M=arma::zeros(n*p,n*p);
  arma::mat B=arma::zeros(n,n);
  arma::mat A=arma::zeros(p,p);
  arma::mat Maux=arma::zeros(n*p,1);
  arma::vec canoni=arma::zeros(n*p);
  
  Progress q(n*p, display_progress);
  for(int k=0;k<n*p;k++) {
    canoni=arma::zeros(n*p);
    canoni(k)=1;
    Maux=arma::zeros(n*p);
    for(int i=0;i<n;i++) {
      
      
      B=K.col(i)*arma::trans(K.col(i));
      A=arma::zeros(p,p);
      for(int j=0;j<numberclose;j++) {
        A=A+wsmall(i,j)*arma::trans(x.row(i)-x.row(neigh(i,j)))*(x.row(i)-x.row(neigh(i,j)));
      }
      // M=M+arma::kron(A,B);
      
      Maux=Maux+armakron(canoni,B,A);
    }
    
    M.col(k)=Maux;
    q.increment(); 
  }
  return M*(2./(n*(n-1)));
}