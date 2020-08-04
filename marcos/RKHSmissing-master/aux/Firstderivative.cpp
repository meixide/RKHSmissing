#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>



// This code was written by C.Meixide-Garcia

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
    
//Compute first derivative of objective function
    arma::mat Firstderivative(int m,int q,int numberclose,arma::mat X,arma::vec Y,arma::mat wsmall,arma::mat K,arma::mat alpha,arma::mat M,int lr,
                              arma::umat neigh) {
      int n=X.n_rows;
      float n1=n;
      
      int l=arma::as_scalar(lr)-1;
      float divi=1./(n1*(n1-1));
      
      arma::mat ycol=arma::repmat(Y,1,n);
      arma::mat xcol=arma::repmat(X.col(l),1,n);
      
     arma::mat coef=arma::zeros(n,n);
    
     arma::mat aux=(2*divi*(ycol-arma::trans(ycol))%(xcol-arma::trans(xcol)));
      for(int i=0;i<n;i++) {
       for(int j=0;j<numberclose;j++) {
         coef(i,neigh(i,j))=wsmall(i,j)*aux(i,neigh(i,j));
        }
      }
      
      
      
      
      arma::mat U=repmat(trans(K.row(0)),1,n)*trans(coef.row(0));
      for (int i=1;i<n;i++){
        U=U+(repmat(trans(K.row(i)),1,n)*trans(coef.row(i)));
      }
      
      int nr=arma::as_scalar(X.n_rows);
      int rmin=((l)*nr);
      int rmax=((l+1)*nr-1);
       
      
      return M.rows(rmin,rmax)*alpha-U;
    }


    