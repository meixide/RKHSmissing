// This code was written by C.Meixide-Garcia

#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
/*vec coef(mat X,vec Y,mat W,double lambda,vec mu,double sigma) {
  
  int i=0;
  int n=X.n_rows;
  int p=X.n_cols;
  vec alphatilde=zeros(n);
  vec beta=zeros(p);*/
  
Rcpp::List wildHSIC(mat X,mat Y,int TestType=1,float Alpha=0.05, Kernel_X,Kernel_Y,int NumBootstrap=300)
    addpath('util') {
    
// okargs =   {'TestType','Alpha', 'Kernel_X','Kernel_Y' ,'WildBootstrap','NumBootstrap'};
 // defaults = {true,0.05, rbf_dot(X),rbf_dot(Y),@bootstrap_series_2, 300};
  //[test_type,alpha, kernel_X,kernel_Y,wild_bootstrap, numBootstrap] = ...
   // internal.stats.parseArgs(okargs, defaults, varargin{:});  
  
  m=X.n_rows();
  //assert(m == Y.n_rows());
 //   assert(test_type==1 || test_type==2)
    
    K = kernel_X(X,X);
  L = kernel_Y(Y,Y);
  H = eye(m)-1/m*ones(m,m);
  Kc = H*K*H;
  Lc = H*L*H;
  statMatrix = Kc.*Lc;
  
  results = bootstrap_null(m,numBootstrap,statMatrix,alpha,wild_bootstrap,test_type);
  // processes = wild_bootstrap(m,numBootstrap);
  
  int ln=20;
  float ar = exp(-1/ln);
  float variance = 1-exp(-2/ln);
  
  mat  A = randn<mat>(length,numBootstrap);
  
  float w=sqrt(variance)*A;
  a = [1 -ar];
  //processes=filter(1,a,w);
  
  
  testStat = m*mean2(statMatrix);
  
  testStats = zeros(numBootstrap,1);
  for process = 1:numBootstrap
    mn = mean(processes(:,process));
  if test_type==1
  matFix = (processes(:,process)-mn)*(processes(:,process)-mn)';
  else
    matFix = processes(:,process)*processes(:,process)';
  end
    testStats(process) =  m*mean2(statMatrix.*matFix );
  end
    
    results.testStat = testStat;
  results.quantile = quantile(testStats,1-alpha);
  results.reject = testStat > results.quantile;
  end
    
  end
    }
    
    
    

    
   
   //Radial basis function inner product
   //Arthur Gretton
   
   //Pattern input format : [pattern1 ; pattern2 ; ...]
   //Output : p11*p21 p11*p22 ... ; p12*p21 ...
   
   mat rbf_dot(mat patterns) {
     sigma = median_heur([patterns]);
     /*if isnan(sigma)
      sigma = 1;
      warning('median heuristic failed')
      end*/
     // kernel = rbf_dot_deg(X,Y,sigma);
     
     //function [H]=rbf_dot_deg(X,Y,deg)
     
     size1=size(X);
     size2=size(Y);
     
     G = sum((X%X),2);
     H = sum((Y%Y),2);
     
     Q = repmat(G,1,Y._nrow());
     R = repmat(H.t(),X.n_row(),1);
     
     H = Q + R - 2*X*Y.t();
     
     
     H=exp(-H/2/deg^2);
     
     return(H)
       
   }
  

      
      