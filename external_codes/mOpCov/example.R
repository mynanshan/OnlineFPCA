library(MASS)

PHI1=function(v,q){
  re=NULL
  for (k in 1:length(v)) {
    re=cbind(re,sqrt(2)*cos(pi*v[k]*(1:q)))
  }
  return(re)
}

####### generate a two dimensional demo data
generate.data = function(n,m,sigma,q=c(3,2)){
  S1=NULL
  S2=NULL
  evl= 1/(1:(prod(q)))^2
  fullCOE=mvrnorm(n,mu=rep(0,prod(q)),Sigma = diag(evl)) 
  Ucoe = order(kronecker(1:q[1],1:q[2]))
  evl = evl[Ucoe]
  
  Z=rep(0,n*m)
  subject = rep(0, n*m)
  for (i in 1:n) {
    s1=runif(m)
    s2=runif(m)
    S1=c(S1,s1)
    S2=c(S2,s2)
    coe=khatri_rao(PHI1(s1,q[1]),PHI1(s2,q[2]))
    z=colSums(apply(coe,2,function(x) x*fullCOE[i,]))+rnorm(m,0,sigma)
    Z[((i-1)*m+1):(i*m)]=z
    subject[((i-1)*m+1):(i*m)] = rep(i,m)
  }
  return(list(location = cbind(S1,S2),Z=Z, subject = subject, q =q)) 
}

