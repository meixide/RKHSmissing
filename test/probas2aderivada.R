
d=distance_mf(X)
wfull=full_Findweight(n,p,X,Y,res$tau,numberclose,wcens)

tau=11 #por probar
x=X
y=Y
tau=res$tau
K=KernelG(x,tau)
neigh=findneigh(n,numberclose,d)
wsmall=findwneigh(n,numberclose,wfull,neigh) 
Secondderivative(n,p,numberclose,x,y,wsmall,K,neigh)

wsmall=findwneigh(n,numberclose,wfull,res$neigh)
alpha=matrix(runif(n*p),ncol=1)
Firstderivative(n,p,numberclose,X,Y,wsmall,res$K,alpha,res$M,1,res$neigh)

