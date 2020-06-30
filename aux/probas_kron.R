n=10
p=5
t=proc.time()
grande(n*numberclose,diag(p)+0.75,diag(n)+0.69)
tempo=proc.time()-t

A=matrix(runif(p*p),nrow=p)
B=matrix(runif(n*n),nrow=n)
t=proc.time()

tempo=proc.time()-t

kronecker(A,B)-kron_for(A,B)



