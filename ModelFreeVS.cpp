// This code was written by C.Meixide-Garcia

#include <Rcpparmadillo.h>
#include <cmath>
#include <stdio.h>

using namespace Rcpp;


//para cout
#include <iostream> 
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp> //para monitorizar o tempo
#include <progress_bar.hpp>

//timer
#include <chrono> 
using namespace std::chrono; 

// Use auto keyword to avoid typing long 
// type definitions to get the timepoint 
// at this instant use function now() 




//Computes a matrix formed by euclidean distance between elements of the input vector
//This is the previous step to calculate the gaussian kernel 
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

//Calculate wij matrix (see paper for more details)

// [[Rcpp::export]]
arma::mat Findweight(int n,int p,arma::mat X,arma::vec Y,float tau,int numberclose,arma::mat wcens) {
  
  /*arma::vec Y,*/
  arma::mat w(n,n);
  arma::mat weight(n,n);
  weight.zeros();
  arma::rowvec d1(n);
  w.zeros();
  weight=w;
  int i=0;
  float cutoff=0;
  float dij=0;
  int j=0;
  float invtau2=1./pow(tau,2);
  
  arma::mat d=arma::zeros(n,n);
  d=arma::pow(distance_mf(X),2.);
  for (i=0;i<n;i++) {
    d(i,i)=arma::datum::inf;
    
  }
  for (i=0;i<n;i++) {
    d1=d.row(i);
    d1=arma::sort(d1);
    cutoff=d1(numberclose-1);
    for (j=0;j<n;j++) {
      dij=d(i,j);
      if (dij<=cutoff) {
        
        weight(i,j)=exp(-dij*invtau2);
      }
    }
  }
  w=weight;
  w=weight%wcens;
  
  
  
  return(w);
  
}

//Compute gaussian kernel 
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



//Compute coefficients for penalization
//Compute largest eigenvalue (see details about optimization on the paper)
arma::mat largesteigenvalue(arma::mat M,int n,int p) {
  
  
  arma::mat A(n,n,arma::fill::zeros);
  arma::mat R(n,M.n_rows);
  arma::vec autoval(n);
  arma::vec gamma(p);
  
  
  for (int i=1;i<p+1;i++) {
    R=M.rows(((i-1)*n),(i*n-1));
    A=R.cols(((i-1)*n),(i*n-1));
    autoval=arma::eig_sym(A);
    gamma(i-1)=max(autoval);
  }
  return(gamma);
  
}

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





arma::vec RKHS(int n,int p,int numberclose,arma::mat X,arma::vec Y,arma::mat w,arma::mat K,arma::mat M,arma::vec tuning,int max,arma::umat neigh) {
  
  
  
  arma::vec alpha=arma::zeros(n*p);
  arma::vec gamma=arma::zeros(p);
  arma::vec b1=arma::zeros(n*p);
  arma::vec b2=arma::zeros(n*p);
  arma::vec alphal=arma::zeros(n);
  arma::vec delta=arma::zeros(n);
  arma::vec fst=arma::zeros(n);
  arma::vec partone=arma::zeros(n);
  arma::vec alphanew=arma::zeros(n);
  float parttwo=0;
  float tun=0;
  int l=1;
  
  gamma=largesteigenvalue(M,n,p);
  
  
  double eps1=10;
  int iteration=0;

  while (eps1>pow(10,-6) && iteration<=max) {
    b1=alpha;
    
    for (l=1;l<p+1;l++) {
      
      alphal=alpha.subvec(((l-1)*n),(l*n)-1);
      fst=Firstderivative(n,p,numberclose,X,Y,w,K,alpha,M,l,neigh);
      delta=alphal*gamma(l-1)-fst;
      partone=delta/gamma(l-1);
      tun=tuning(l-1);
      parttwo=1-tun/pow(dot(delta,delta),(1./2));
      if(parttwo<0) {
        parttwo=0;
        
      }
      
      alphanew=partone*parttwo;
      alpha.subvec(((l-1)*n),(l*n-1))=alphanew;
      
    }
    
    b2=alpha;
    
    eps1=pow(dot((b1-b2),(b1-b2)),(1./2));
    iteration=iteration+1;
    
    
  }
  return alpha; //hai que devolver alpha
  
  
}





//Compute stability of varaible selection method using cross validation 
arma::mat Crossvalidation(arma::mat X,arma::vec Y,int n,int p,int n1,int n2,int numberclose,int totaltime,arma::vec alphatildeKnorm,arma::mat wcens,
                          bool display_progress=true) {
  //It is used to choose the tuning parameter based on training data
  
  float t=20;
  int i=0;
  int j=0;
  int t0=0;
  arma::vec alpha1(n1*p);
  arma::vec alpha2(n2*p);
  //The dimension choose for each data set
  arma::colvec criteria1(p,arma::fill::zeros);
  arma::colvec criteria2(p,arma::fill::zeros); 
  arma::mat STAB(totaltime,t,arma::fill::zeros);
  arma::cube CRITERIA1(p,t,totaltime,arma::fill::zeros);
  arma::cube CRITERIA2(p,t,totaltime,arma::fill::zeros);
  arma::mat slicek(p,t);
  arma::mat Distance1(n1,n1,arma::fill::zeros);
  arma::uvec xulie(n);
  arma::mat K1(n1,n1,arma::fill::zeros);
  arma::mat M1(n1*p,n1*p);
  arma::mat w1(n1,n1);
  arma::mat w1full(n1,n1);
  arma::mat w2(n2,n2);
  arma::mat w2full(n2,n2);
  arma::mat M2(n2*p,n2*p);
  float lambda;
  
  //The step length and start point for tuning parameter
  int start=-2;
  int step=4;
  //Do the cross validation for totaltime times
Progress q(totaltime, display_progress);
  for (int k=0;k<totaltime;k++) {
    //Randomly divide the data into two parts
    arma::mat X1(n1,p,arma::fill::zeros);
    arma::mat Y1(n1,1,arma::fill::zeros);
    arma::mat X2(n2,p,arma::fill::zeros);
    arma::mat Y2(n2,1,arma::fill::zeros);
    
    // set some values:
    
    for (int i=0; i<n; i++) xulie(i)=i; // 0 1 2 3 4 5 6 7 8 9
    
    // using built-in random generator:
    std::random_shuffle ( xulie.begin(), xulie.end() );
    
    
    for (j=0;j<n1;j++) {
      X1.row(j)=X.row(xulie(j));
      Y1.row(j)=Y.row(xulie(j));
    }
    
    for (j=0;j<n2;j++) {
      X2.row(j)=X.row(xulie(j+n1));
      Y2.row(j)=Y.row(xulie(j+n1));
    }
    
    
    //Find the Gaussian Kernel matrix
    
    
    Distance1=distance_mf(X1);
    
    arma::vec med(n1,arma::fill::zeros);
    for(i=0;i<n1;i++){
      med(i)=arma::median(Distance1.col(i));
    }
    float tau1=arma::median(med);
    K1=KernelG(X1,tau1);
    arma::mat Distance2(n2,n2,arma::fill::zeros);
    arma::mat K2(n2,n2,arma::fill::zeros);
    Distance2=distance_mf(X2);
    arma::vec med2(n2,arma::fill::zeros);
    for(i=0;i<n2;i++){
      med2(i)=arma::median(Distance2.col(i));
    }
    float tau2=arma::median(med2);
    K2=KernelG(X2,tau2);
    
    //Compute the weight and decide how many close points will be selected for these two data sets
    
    
    
    arma::mat wcens1(n1,n1);
    arma::mat wcens2(n2,n2);
    arma::uvec rowscols1(n1);
    arma::uvec rowscols2(n2);
    rowscols1=xulie.subvec(0,(n1-1));
    rowscols2=xulie.subvec(n2,(2*n2-1));
    wcens1=wcens.submat(rowscols1,rowscols1);
    wcens2=wcens.submat(rowscols2,rowscols2);
    
   
    arma::mat wfull1=full_Findweight(n1,p,X1,Y1,tau1,numberclose,wcens1);
    arma::umat neigh1=findneigh(n1,numberclose,Distance1);
    arma::mat wsmall1=findwneigh(n1,numberclose,wfull1,neigh1);
    M1=Secondderivative(n1,p,numberclose,X1,Y1,wsmall1,K1,neigh1);  
    
    arma::mat wfull2=full_Findweight(n2,p,X2,Y2,tau2,numberclose,wcens2);
    arma::umat neigh2=findneigh(n2,numberclose,Distance2);
    arma::mat wsmall2=findwneigh(n2,numberclose,wfull2,neigh2);
    M2=Secondderivative(n2,p,numberclose,X2,Y2,wsmall2,K2,neigh2);
    
    int max=10;
    arma::colvec tuning1(p,arma::fill::zeros);
    arma::colvec tuning2(p,arma::fill::zeros);
    
    for (t0=0;t0<t;t0++) {
      arma::colvec VS1(p,arma::fill::zeros);
      arma::colvec VS2(p,arma::fill::zeros);
      lambda=pow(10.,start+(t0)*step/19.);
      //The adaptive tuning parameter
      for (i=0;i<p;i++) {
        tuning1(i)=lambda/alphatildeKnorm(i);
        tuning2(i)=lambda/alphatildeKnorm(i);
      }
      //Get the estimation for alpha for these two data sets
      alpha1=RKHS(n1,p,numberclose,X1,Y1,wsmall1,K1,M1,tuning1,max,neigh1);
      alpha2=RKHS(n2,p,numberclose,X2,Y2,wsmall2,K2,M2,tuning2,max,neigh2);
      
      for (i=0;i<p;i++) {
        
        
        criteria1(i)= arma::norm(alpha1.subvec((i)*n1,((i+1)*n1)-1),2);
        criteria2(i)=arma::norm(alpha2.subvec((i)*n2,((i+1)*n2)-1),2);
        if (criteria1(i)>=0.001) {
          VS1(i)=1;
        }
        
        if (criteria2(i)>=0.001) {
          VS2(i)=1;
        }
        
      }
      int n11=0;
      int n12=0;
      int  n21=0;
      int n22=0;
      for (i=0;i<p;i++) {
        if (VS1(i)==1 && VS2(i)==1) {
          n11=n11+1;
        }
        
        if (VS1(i)==1 && VS2(i)==0) {
          n12=n12+1;
        }
        
        if (VS1(i)==0 && VS2(i)==1) {
          n21=n21+1;
        }
        
        if (VS1(i)==0 && VS2(i)==0) {
          n22=n22+1;
        }
      }
      //Find the variable selection stability based on these two data sets
      float Pra=(n11+n22)/(float)p;
      float Pre=(float)(n11+n12)*(n11+n21)/(pow(p,2))+(float)(n12+n22)*(n21+n22)/(pow(p,2));
      STAB(k,t0)=(Pra-Pre)/(1-Pre);
      if (sum(VS1)==p) {
        STAB(k,t0)=-1;
      }
      
      if (sum(VS2)==p) {
        STAB(k,t0)=-1 ;
      }
      
      if (sum(VS1)==0) {
        STAB(k,t0)=-1;
      }
      if (sum(VS2)==0) {
        STAB(k,t0)=-1 ;
      }
      
      CRITERIA1.subcube(0,t0,k,p-1,t0,k)=criteria1;
      
      CRITERIA2.subcube(0,t0,k,p-1,t0,k)=criteria2;    }
    q.increment(); 
  
  }
  return STAB;
}

//Compute the best choice of variables to explain the model
arma::vec findDIM(arma::mat X,arma::vec Y,int n,int p,float lambda,int numberclose,arma::vec alphatildeKnorm,arma::mat wcens) {
  //It is used to find the useful dimensions
  arma::mat K(n,n,arma::fill::zeros);
  arma::vec alpha=arma::zeros(n*p);
  
  arma::mat w(n,n,arma::fill::zeros);
  arma::mat wfull(n,n,arma::fill::zeros);
  arma::mat M(n*p,n*p,arma::fill::zeros);
  //Find the Gaussian Kernel matrix
  arma::mat Distance(n,n,arma::fill::zeros);
  int i=0;
  Distance=distance_mf(X);
  
  arma::mat wcens3(n,n);
  wcens3=wcens.submat(0,0,n-1,n-1);
  
  arma::vec med(n,arma::fill::zeros);
  for(i=0;i<n;i++){
    med(i)=arma::median(Distance.col(i));
  }
  float tau=arma::median(med);
  
  
  K=KernelG(X,tau);
  //Find the weight for each training data point
  wfull=full_Findweight(n,p,X,Y,tau,numberclose,wcens3);
  arma::umat neigh=findneigh(n,numberclose,Distance);
  arma::mat wsmall=findwneigh(n,numberclose,wfull,neigh);
  M=Secondderivative(n,p,numberclose,X,Y,wsmall,K,neigh);
 
  
  arma::colvec tuning(p,arma::fill::zeros);
  
  for (i=0;i<p;i++) {
    tuning(i)=lambda/alphatildeKnorm(i);
  }
  int max=50;
  //Find the unknown parameter alpha based on given tuning parameter
  alpha=RKHS(n,p,numberclose,X,Y,wsmall,K,M,tuning,max,neigh);
  //Find the norm for each part of alpha
  arma::colvec criteria(p,arma::fill::zeros);
  for (i=0;i<p;i++) {
    criteria(i)=pow(arma::dot(alpha.subvec(i*n,(i+1)*n-1),alpha.subvec(i*n,(i+1)*n-1)),(1./2));
  }
  return criteria;
}

// [[Rcpp::export]]

//Main program
Rcpp::List ModelFreeVS(arma::mat X,arma::vec Y,int numberclose,int totaltime,arma::mat wcens) {
  //------------------------------------------------------------------
  //------------------------------------------------------------------
  //------------------------------------------------------------------
  //------------------------------------------------------------------
  //Find the alphatilde to give the adaptive lasso weight
  
  int p=X.n_cols;
  arma::vec criteria(p);
  int n=X.n_rows;
  
  int n1=0;
  int n2=0;
  arma::mat K(n,n,arma::fill::zeros);
  arma::mat w(n,n,arma::fill::zeros);
  arma::mat wfull(n,n,arma::fill::zeros);
  
  arma::mat M(n*p,n*p,arma::fill::zeros);
  arma::Col<double> betatilde=arma::zeros(n*p,1);
  
  int i=0;
  
  
  //Find the Gaussian Kernel matrix

  arma::mat Distance(n,n,arma::fill::zeros);
  Distance=distance_mf(X);

 
  
  // Subtract stop and start timepoints and 
  // cast it to required unit. Predefined units 
  // are nanoseconds, microseconds, milliseconds, 
  // seconds, minutes, hours. Use duration_cast() 
  // function. 
 
 //***********************************************************
  
  auto start = high_resolution_clock::now(); 
  arma::vec med(n,arma::fill::zeros);
  for(i=0;i<n;i++){
    med(i)=arma::median(Distance.col(i));
  }
  float tau=arma::median(med);
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<milliseconds>(stop - start); 
  
  Rcout << "Time taken by function median: "
        << duration.count() << " millseconds" << endl; 
  
  //***********************************************************
  K=KernelG(X,tau);
  auto startff = high_resolution_clock::now(); 
  
  wfull=full_Findweight(n,p,X,Y,tau,numberclose,wcens);
  
  auto stopff = high_resolution_clock::now(); 
  auto durationff = duration_cast<milliseconds>(stopff - startff); 
  
  Rcout << "Time taken by function full_Findweight: "
        << durationff.count() << " millseconds" << endl; 
  
  //***********************************************************
  
  
  arma::umat neigh=findneigh(n,numberclose,Distance);
 arma::mat wsmall=findwneigh(n,numberclose,wfull,neigh);
 
 //w=Findweight(n,p,X,Y,tau,numberclose,wcens);

  
  //***********************************************************
  
 

  auto startsd = high_resolution_clock::now(); 

  
  M=Secondderivative(n,p,numberclose,X,Y,wsmall,K,neigh);
  auto stopsd = high_resolution_clock::now(); 
  auto durationsd = duration_cast<milliseconds>(stopsd - startsd); 
  
  Rcout << "Time taken by function Secondderivative: "
        << durationsd.count() << " millseconds" << endl; 
  

  arma::Col<double> tuning(p,arma::fill::zeros);
  tuning=tuning+0.000001;
  auto startr = high_resolution_clock::now(); 
  
  //50 iteration to get estimation for alphatilde
  betatilde=RKHS(n,p,numberclose,X,Y,wsmall,K,M,tuning,50,neigh);
  
  auto stopr = high_resolution_clock::now(); 
  auto durationr = duration_cast<milliseconds>(stopr - startr); 
  
  Rcout << "Time taken by function RKHS: "
       << durationr.count() << " milliseconds" << endl; 
  
  
  
  arma::colvec alphatildeKnorm(p,arma::fill::zeros);
  float atkn=0;
  
  for (i=0;i<p;i++) {
    atkn=dot(betatilde.subvec(i*n,(i+1)*n-1),K*betatilde.subvec(i*n,(i+1)*n-1));
    alphatildeKnorm(i)=atkn;
  }
  //*********ata aqui todo ben 
  
  //Variable selection stability for a sequence of tuning parameters
  
  n1=floor(n/2.);
  n2=n1;
  
  
  
  arma::mat STAB=Crossvalidation(X,Y,n,p,n1,n2,numberclose,totaltime,alphatildeKnorm,wcens);
  float t=0.;
  
  
  arma::vec biaozhun(20,arma::fill::ones);
  
  arma::vec ms(20);
  ms=arma::trans(arma::mean(STAB,0));
  
  if (arma::max(ms)>0) {
    biaozhun=biaozhun*arma::max(ms)*0.9;
    arma::uvec tt1=find(ms>=biaozhun);
    t=tt1(0);
    
  }
  
  if (arma::max(ms)<=0) {
    biaozhun=biaozhun*arma::max(ms);
    arma::uvec tt2=find(ms==biaozhun);
    t=tt2(0);
    
  }
  
  t=t+1; 
  
  
  //Decide which variable is informative
  float lambda=pow(10.,(-2.+(t-1)*4./19));
  arma::mat X3=X.submat(0,0,n1-1,p-1);
  arma::vec Y3=Y.subvec(0,n1-1);
  
  
  criteria=findDIM(X3,Y3,n1,p,lambda,numberclose,alphatildeKnorm,wcens);
  
  
  
   Rcpp::List ret;
   ret["K"] = K;
   ret["M"] = M;
   ret["w"] = w;
   ret["floor"] = floor(n/2.);
   ret["alphatildeKnorm"] = alphatildeKnorm;
   ret["STAB"] = STAB;
   ret["t"] = t;
   ret["tau"]=tau;
   ret["criteria"] = criteria;
   ret["Distance"]=Distance;
   ret["X3"] = X3;
   ret["Y3"] = Y3;
   ret["ms"] = ms;
   ret["lambda"] = lambda;
   ret["betatilde"] = betatilde;
   ret["neigh"] = neigh;
   ret["wsmall"] = neigh;
   
   
  
  return ret;
}