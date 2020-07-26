library("SuperLearner")
library(copula)
library(psych)
library(gtools)
library("ranger")
library("glmnet")
library("gam")
library("kernlab")

n=200
myCop <- normalCopula(param=c(0.3,0.7,0.8), dim = 3, dispstr = "un")
myMvd <- mvdc(copula=myCop, margins=c("unif", "unif", "unif"),
              paramMargins=list(list(min=0,max=10),
                                list(min=0,max=10), 
                                list(min=0,max=10)) )


Z <- rMvdc(n, myMvd)

eps=rnorm(n,0,1)
Y=2*(Z[,1])+eps
colnames(Z) <- c("x1", "x2", "x3")
pairs.panels(Z)


beta=c(-0.05,0,0)
p=numeric(n) #vector with pi
for(i in 1:n) {
  x=Z[i,]
  p[i]=inv.logit(sum(x*beta), min = 0, max = 1)
}


M=numeric(n)
for(i in 1:n) {
  M[i]=rbinom(1,1,1-p[i])
}
Yorig=Y
Y[!M]=NA


res=incomplete.reg(Z,Y)


plot(p~as.numeric(res$p))
plot(Y[!is.na(Y)],res$Ypred[!is.na(Y)])


