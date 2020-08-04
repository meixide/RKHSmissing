Y[is.na(Y)]=0
mu=res$Ypred
lambda=0.001

d=as.numeric(M/p)

Ytilde=Y[res$p>0.1]
Ztilde=Z[res$p>0.1,]
mutilde=mu[res$p>0.1]
dtilde=d[res$p>0.1]


#d=d/sum(d)
Wtilde=(diag(dtilde))
W=diag(d)
#alphatilde=coef(Ztilde,Ytilde,Wtilde,lambda,mutilde,3)
alphatilde=coef(Z,Y,W,lambda,mu,15)
prediction=numeric(n);

for(i in 1:n) {
  prediction[i]=kernel_machine(Z,as.numeric(Z[i,]),as.numeric(alphatilde),15)
}

plot(Yorig~mu)
plot(Yorig,prediction)
