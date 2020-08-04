library("SuperLearner")

n=200
p=5

x1=runif(n,0,10)
x2=runif(n,0,10)
x3=runif(n,0,10)
x4=runif(n,0,10)
x5=runif(n,0,10)

eps=rnorm(n,0,1)

Y=2*(x1+x2+x3)+(x4+x5)+eps


mod.lm=lm(Y~x1+x2+x3+x4+x5)
summary(mod.lm)

x_train=data.frame(x1,x2,x3,x4,x5)
y_train=Y

library(psych)
pairs.panels(x_train)

library("ranger")
library("glmnet")
sl = SuperLearner(Y = y_train, X = x_train, family = gaussian(),
                  SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"))
sl

library("gam")
library("kernlab")
cv_sl = CV.SuperLearner(Y = y_train, X = x_train, family = gaussian(),
                  
                        V = 5,SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.gam", "SL.ksvm"))

#knn so vale para binomial()

summary(cv_sl)
plot(cv_sl) 




library(copula)

myCop <- normalCopula(param=c(0.3,0.7,0.8), dim = 3, dispstr = "un")
myMvd <- mvdc(copula=myCop, margins=c("unif", "unif", "unif"),
              paramMargins=list(list(min=0,max=10),
                                list(min=0,max=10), 
                                list(min=0,max=10)) )


Z2 <- rMvdc(200, myMvd)
Y=2*(Z2[,1]+Z2[,2])+Z2[,3]
colnames(Z2) <- c("x1", "x2", "x3")
pairs.panels(Z2)

cv_sl = CV.SuperLearner(Y = Y, X = as.data.frame(Z2), family = gaussian(),
                        
                        V = 6,SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.gam", "SL.ksvm"))

summary(cv_sl)
plot(cv_sl) 

sl = SuperLearner(Y = y_train, X = x_train, family = gaussian()
                        
                        ,SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.gam", "SL.ksvm"))

pred = predict(sl, x_train, onlySL = TRUE)

plot(y_train~as.numeric(pred[[1]]))

