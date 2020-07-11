

#executar dende aqui


library("Rcpp")

sourceCpp("ModelFreeVS.cpp")
#install.packages("RcppProgress")
library("RcppProgress") #hai que engadir tres linhas por rutina
setwd("~/repodir/RKHSmissing")

n=200
p=8

pesos=rep(1,n)
wcens= pesos%*%t(pesos);
wcens= wcens*n*(n-1)/2; 

X<-matrix(0,nrow=n,ncol=p)

# for (i in 1:p) {
#   X[,i]<-runif(n,0,10)
# }
# 
# Y<-2*X[,2]


totaltime<-5
numberclose<-4


X=read.csv('datn200p8.dat')
X=as.matrix(X[,-1])

Y<-rep(0,n)
Y<-2*X[,3]


dim(wcens)
dim(X)
t=proc.time()
res<-ModelFreeVS(X,Y,numberclose,totaltime,wcens)
tempo=proc.time()-t
#write.csv(X,'datn200p8.dat')
# library("profvis")
# 
# profvis(
#   res<-ModelFreeVS(X,Y,numberclose,10,wcens)
# ) non distingue entre subrutinas
