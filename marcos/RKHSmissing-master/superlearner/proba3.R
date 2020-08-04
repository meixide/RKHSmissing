library("SuperLearner")
library(copula)
library(psych)
library(gtools)
library("ranger")
library("glmnet")
library("gam")
library("kernlab")
library(WeightIt)



n=2000
p=5

myCop <- normalCopula(rep(0.4,p * (p - 1) / 2), dim = 5, dispstr = "un")
myMvd <- mvdc(copula=myCop, margins=c("unif", "unif", "unif","unif","unif"),
              paramMargins=list(list(min=0,max=10),
                                list(min=0,max=10), 
                                list(min=0,max=10),
                                list(min=0,max=10), 
                                list(min=0,max=10)))


Z <- rMvdc(n, myMvd)
#Z= Z[,1]
eps=rnorm(n,0,1)
#Y= 2*Z[,1]+1.5*Z[,2] + 1*Z[,3] + 0.5*Z[,4]+eps
Y=2*(Z[,1])+eps
#colnames(Z) <- c("x1", "x2", "x3")
#pairs.panels(Z)


beta=c(-0.01,-0.01,-0.01,-0.01,-0.01)
#beta= beta[1]
p=numeric(n) #vector with pi
for(i in 1:n) {
  x=Z[i,]
  #x= Z[i]
  
  p[i]=inv.logit(sum(x*beta), min = 0, max = 1)
}


M=numeric(n)
for(i in 1:n) {
  M[i]=rbinom(1,1,1-p[i])
}
Yorig=Y
Y[!M]=NA
#setwd("~/Descargas/RKHSmissing-master/superlearner")
source("incompletenew.R")


res=incomplete.reg(Z,Y)

plot(p~as.numeric(res$p))
points(p[M==1]~as.numeric(res$p[M==1]),col='red')
plot(Y[!is.na(Y)],res$Ypred[!is.na(Y)])

plot(Yorig,res$Ypred)

####################################################33
#setwd("~/Descargas/RKHSmissing-master/kernel_machine")
library("Rcpp")
library(RcppArmadillo)
sourceCpp("../kernel_machine/regression.cpp")
Y[is.na(Y)]=0
mu=res$Ypred
lambda=0.05




d=as.numeric(M/p)

Ytilde=Y[res$p>0.1]
#Ztilde=Z[res$p>0.1,]

Ztilde=Z[res$p>0.1]

mutilde=mu[res$p>0.1]
dtilde=d[res$p>0.1]


#d=d/sum(d)
Wtilde=(diag(dtilde))
W=diag(d)
#alphatilde=coef(Ztilde,Ytilde,Wtilde,lambda,mutilde,3)
Z= as.matrix(Z)



distanciasZ= dist(Z)
median(distanciasZ^2)

waux= M/res$p
#waux= M/p
w= diag(as.numeric(waux),n,n) 
# borre W
alphatilde=coef(Z,Y,w,0.5,mu,sqrt(1/70.26))
prediction=numeric(n);

for(i in 1:n) {
  prediction[i]=kernel_machine(Z,as.numeric(Z[i,]),as.numeric(alphatilde),sqrt(1/70.26))
}

K=Kernelfull(Z,sqrt(1/70.26)) 
pred2=K%*%alphatilde
plot(pred2,Yorig)
points(pred2[M==1],Yorig[M==1],col='red')
plot(prediction, Yorig)
#plot(alphatilde,new)
# library("listdtr")
# kernelridgestandart=krr(Z, Y, group = NULL)

#library("krr")

lambda= seq(0.5,3,by=0.1)
sigma= seq(0.5,30,by=1)

combinaciones=expand.grid(lambda,sigma)


errores= numeric(dim(combinaciones)[1])
for(i in 1:(dim(combinaciones)[1])){
  aux= krr(Z,Yorig ,lambda = combinaciones[i,1],sigma = combinaciones[i,2])
  errores[i]= aux$MSE
  
}


waux= M/res$p
#waux= M/p
w= diag(as.numeric(waux),n,n) 

rbf <- rbfdot(sigma = 8.5)
K <- kernelMatrix(rbf, Z)
lambda= 0.5
new=solve(w%*%K+diag(lambda,n,n))%*%(w%*%Y)
new2= K%*%new
plot(Yorig,new2)


max(Kernelfull(Z, sqrt(1/8.5))-K)




kernlab::kernelMult(Z,Z)


plot(Yorig~mu)
plot(Yorig,prediction)
