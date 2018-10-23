
#Normalize density function in such a way that max(rho)<=1
exit_density=function(x,k=0.1,r=1){
  return(exp(k*r*(cos(x)-1)))
}



rdens=function(n,rho,xmin=0,xmax=2*pi,params=list()){

  N=5*n
  y=c()
  while(length(y)<n){
  lambda=runif(N)
  x=runif(N,xmin,xmax)
  r=do.call(rho,c(list(x=x),params))
  y=c(y,x[lambda<rho(x)])
  }

  return(y[1:n])
}
