#executar dende aqui


library("Rcpp")
#setwd("C:/Users/Carlos M Garcia/Desktop/mfmissing")


#install.packages("RcppProgress")
library("RcppProgress") #hai que engadir tres li?as por rutina

#setwd("C:/Users/Carlos M Garcia/Google Drive/Marcos/missingdata/R/modelfreecensoredfinal")
#setwd("/home/carlos/R")
n=200
p=8

pesos=rep(1,n)
wcens= pesos%*%t(pesos);
wcens= wcens*n*(n-1)/2; 

X<-matrix(0,nrow=n,ncol=p)




for (i in 1:p) {
 X[,i]<-runif(n,0,10)
 }
# 
# Y<-2*X[,2]


totaltime<-5
numberclose<-4


X=read.csv('datn200p8.dat')
X=as.matrix(X[,-1])

Y<-rep(0,n)
Y<-2*X[,3]


totaltime<-5
numberclose<-4
tau=11.4 #por probar
d=distance_mf(X)
wfull=full_Findweight(n,p,X,Y,tau,numberclose,wcens)


x=X
y=Y

K=KernelG(x,tau)
neigh=findneigh(n,numberclose,d)
wsmall=findwneigh(n,numberclose,wfull,neigh) 
t1=proc.time()
res<-ModelFreeVS(X,Y,numberclose,totaltime,wcens)
t2=proc.time()
M2=Secondderivativevello(n,p,numberclose,x,y,wsmall,K,neigh)
tempo2=proc.time()-t2
t1=proc.time()
M1=Secondderivative(n,p,numberclose,x,y,wsmall,K,neigh)
tempo1=proc.time()-t1
t3=proc.time()
M3=Secondderivativeprobas(n,p,numberclose,x,y,wsmall,K,neigh)
tempo3=proc.time()-t3
t4=proc.time()
M4=Secondderivativearma(n,p,numberclose,x,y,wsmall,K,neigh)
tempo4=proc.time()-t4
S=t(X)%*%X
v=runif(n*p,0,1)


testkron=armakron(v,S,K)
A=matrix(c(1,2,3,4),nrow=2)
B=matrix(c(1,1,1,1),nrow=2)
C=matrix(c(5,5,5,5),nrow=2)

kronecker(A,B+C)
kronecker(A,B)+kronecker(A,C)


################################

wsmall=findwneigh(n,numberclose,wfull,res$neigh)
alpha=matrix(runif(n*p),ncol=1)
Firstderivative(n,p,numberclose,X,Y,wsmall,res$K,alpha,res$M,1,res$neigh)

