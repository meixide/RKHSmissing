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
  X= as.data.frame(X)
  
  pX= dim(X)[2]
  
#   if(pX==1){
#     repetir= rep(1,n)
#     X= data.frame(X,repetir)
#     pX=2
#   }
#   
# X= as.data.frame(X)  

# Añadi GLM


  slmis = SuperLearner(Y = M, X = X, family = binomial(),SL.library = c("SL.mean","SL.glm", "SL.ranger", "SL.gam", "SL.ksvm"))
  
  # porque tiñas Z
  
  pred = predict(slmis, as.data.frame(X), onlySL = TRUE)
  p=pred[[1]]
  
  Ycomplete=Y[!is.na(Y)]
  Xnew= X[!is.na(Y),]
  Xnew= as.data.frame(Xnew)
  colnames(Xnew)= colnames(X)
  if(pX>1){
  
  sl = SuperLearner(Y = Ycomplete, X = as.data.frame(Xnew), family = gaussian(),SL.library = c("SL.mean","SL.lm", "SL.ranger", "SL.gam", "SL.ksvm"))
  
  pred2 = predict(sl, newdata=as.data.frame(X))
  Ypred=pred2[[1]]
  
  }else{
    
    sl = SuperLearner(Y = Ycomplete, X = Xnew, family = gaussian(),SL.library = c("SL.mean","SL.ranger", "SL.lm","SL.gam", "SL.ksvm"))
    pred2 = predict(sl, newdata=as.data.frame(X))
    Ypred=pred2[[1]]
    
    
  }
  
  
  return(list("Ypred" = Ypred , "p" = p))
  
}