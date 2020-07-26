library("SuperLearner")
library(copula)
library(psych)
library(gtools)
library("ranger")
library("glmnet")
library("gam")
library("kernlab")

incomplete.reg=function(X,Y) {
  n=length(Y)
  M=as.numeric(is.na(Y))
  
  slmis = SuperLearner(Y = M, X = as.data.frame(X), family = binomial()
                    
                    ,SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.gam", "SL.ksvm"))
  pred = predict(slmis, as.data.frame(Z), onlySL = TRUE)
  p=pred[[1]]
  
  Ycomplete=na.omit(Y)
  
  sl = SuperLearner(Y = Ycomplete, X = as.data.frame(X[!is.na(Y),]), family = gaussian()
                    
                    ,SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.gam", "SL.ksvm"))
  
  pred2 = predict(sl, as.data.frame(X), onlySL = TRUE)
  Ypred=pred2[[1]]
  
  return(list("Ypred" = Ypred , "p" = p))
  
}