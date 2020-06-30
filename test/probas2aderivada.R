
d=distance_mf(X)
wfull=full_Findweight(n,p,X,Y,res$tau,numberclose,wcens)


res$M-Secondderivative(n,p,numberclose,X,Y,w,res$K)

wsmall=findwneigh(n,numberclose,wfull,res$neigh)
alpha=matrix(runif(n*p),ncol=1)
Firstderivative(n,p,numberclose,X,Y,wsmall,res$K,alpha,res$M,1,res$neigh)

