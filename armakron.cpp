#include <Rcpparmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

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