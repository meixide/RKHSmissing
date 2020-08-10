// This code was written by C.Meixide-Garcia

#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]

//Dino Sejdinovic, 2013
//D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.

//median heuristic bandwidth selection

float median_heur(mat Z) {
  
  mat Zmed;
  
  int size1=Z.n_rows;
  if (size1>100){//choose 100 random samples from Z if it contains
  //more than 100 points
  uvec ind=sort_index(randu(size1));
  Zmed = Z.cols(ind(0),ind(99));
  size1 = 100;
  }else{
  Zmed = Z;
    }
  
  mat G = sum((Zmed%Zmed),1);
  mat Q = repmat(G,0,size1);
  mat R = repmat(G.t(),size1,0);
  mat dists = Q + R - 2*(Zmed*Zmed.t());
  dists = dists-trimatl(dists);
  vec distsr = reshape(dists,size1^2,1);
  distsr = sqrt(distsr);
  float sig = sqrt(0.5*median(dists(find(dists>0))));
  return(sig);
  
  
}

