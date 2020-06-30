#executar dende aqui


library("Rcpp")
setwd("C:/Users/Carlos M Garcia/Desktop/mfmissing")

sourceCpp("ModelFreeVS.cpp")
#install.packages("RcppProgress")
library("RcppProgress") #hai que engadir tres li?as por rutina

#setwd("C:/Users/Carlos M Garcia/Google Drive/Marcos/missingdata/R/modelfreecensoredfinal")
#setwd("/home/carlos/R")
n=200
p=8

pesos=rep(1,n)
wcens= pesos%*%t(pesos);
wcens= wcens*n*(n-1)/2; 

# X<-matrix(0,nrow=n,ncol=p)
# 
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

res$criteria

#se todo funciona ben deberia dar isto

# [,1]
# [1,] 0.0000000000
# [2,] 0.0000000000
# [3,] 0.5093863934
# [4,] 0.0000000000
# [5,] 0.0000000000
# [6,] 0.0000000000
# [7,] 0.0000000000
# [8,] 0.0006083079


#write.csv(X,'datn200p8.dat')
# library("profvis")
# 
# profvis(
#   res<-ModelFreeVS(X,Y,numberclose,10,wcens)
# ) non distingue entre subrutinas
