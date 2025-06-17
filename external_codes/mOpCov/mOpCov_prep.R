library(RSpectra)
library(MASS)

getM.sam <- function(location,nsam,ker){
  d = length(q)
  M = matrix(nrow = nrow(location),ncol = sum(q))
  a = c(0,cumsum(q))
  for (i in 1:d) {
    if(ker == "sob"){ K <- getK_sob(location[,i])}
    if(ker == "cos"){ K <- getK_cos(location[,i])}
    ind <- sample(length(location[,i]), nsam)
    u <- svd(K[,ind])$u
    a <- chol(t(u) %*% K %*% u)
    M[,((a[i]+1):a[i+1])] = Eig$vectors%*%diag(sqrt(Eig$values))
  }
  return(M)
}

getM<- function(location, q, tol=1e-6,ker){
  d = length(q)
  M = matrix(nrow = nrow(location),ncol = sum(q))
  a = c(0,cumsum(q))
  for (i in 1:d) {
    if(ker == "sob"){ K <- getK_sob(location[,i])}
    if(ker == "cos"){ K <- getK_cos(location[,i])}
    Eig=eigs_sym(K,q[i])
    M[,((a[i]+1):a[i+1])] = Eig$vectors%*%diag(sqrt(Eig$values))
  }
  return(M)
}

prep <- function(location, x, subject, q, Mmethod="eig", tol=1e-6, nsam=50, ker="sob"){
  d = length(q)

  Xs <- list()
  ms <- NULL
  i <- 0
  for (zz in unique(subject)){
    if (sum(subject==zz)>1){
      i <- i+1
      Xs[[i]] <- as.double(x[subject == zz])
      ms <- c(ms, length(Xs[[i]]))
    }
  }
  n <- length(Xs)
  if (Mmethod=="eig"){
    M <- getM(location = location,q = q,ker=ker)
  }  else if (Mmethod=="sam"){
    M <- getM.sam(location,nsam=nsam,ker)
  }
  ii <- c(0, cumsum(ms))
  Ms <- list()
  for (i in (1:n)){
      Mstemp = M[(ii[i]+1):ii[i+1],]
    Ms[[i]] = Mstemp
  }
  
  obj <- Qlossprep_cpp(Xs, Ms, q, 1:n)
  return(list(Xs=Xs, Ms=Ms, ms=ms, M=M, q=q, location = location, G=obj$G, h=obj$h, c=obj$c))
}


mOpCov.control <- function(Mmethod="eig", nsam=50, eta=1e-07, maxiter=6000, preptol=1e-8,tol=1e-10){
  return(list(Mmethod=Mmethod,  nsam=nsam, preptol=preptol, eta=eta, maxiter=maxiter, tol=tol))
}

##################################
#########  mOpCov ################
##################################

mOpCov = function(location, x, subject,q, lam, b0=NULL, ker = "sob", control=list()){
  control <- do.call(mOpCov.control, control)
  prepout = prep(location = location,x = x,subject =  subject,q = q,Mmethod=control$Mmethod, tol=control$preptol, nsam=control$nsam,ker=ker)
  if (is.null(b0)){b = rep(0, prod(q^2))}else{b=b0}
  
  eig = eigen(prepout$G)
  GU =eig$vectors
  vG = eig$values
  if(control$eta>10*lam$lam) {eta = 10*lam$lam}else{eta = control$eta}
  res = ACADMM(q,control$maxiter,prepout$G,GU,vG,prepout$h,b,control$eta,lam$lam*lam$alpha,lam$lam*(1-lam$alpha)/2,control$tol)
  
  Minv = list()
  a = c(0,cumsum(q))
  for (i in 1:length(q)) {
   Minv[[i]] =  ginv(prepout$M[,((a[i]+1):a[i+1])])
  }

  b = res$b
  d = length(q)
  perm = c(rev(1:d),rev((1+d):(2*d)))
  Bt = as.tensor(array(unfold(b,perm,c(rev(q),rev(q))),c(q,q)))
  r = rep(0,d)
  for (i in 1:d) {
    Mtemp = k_unfold(Bt,m=i)@data
    r[i] = rankMatrix(Mtemp,tol=1e-5)[1]
  }
  
  if(!all(diff(r-q)<0)) {
    for(i in 1:20){
      Bt = as.tensor(array(unfold(b,perm,c(rev(q),rev(q))),c(q,q)))
      Bt_low = tucker(Bt,c(r,r))
      Bm = matrix(unfold(as.vector(Bt_low$est@data),perm,c(q,q)),nrow=prod(q))
      Bm = 0.5 * Bm +0.5 *t(Bm)
      # ensure pd
      eig = eigen(Bm)
      ind = eig$values>=0
      Bm = eig$vectors[,ind] %*% diag(eig$values[ind]) %*% t(eig$vectors[,ind])
      b = as.vector(Bm)
    }
    R = rankMatrix(Bm,tol=1e-10)[1]
    res$b = b
  }
  return(list(solution=res$b,Minv=Minv,q=q,location=prepout$location,ker = ker, R= R, r =r))
}


mOpCov.cv.control <- function(lams=NULL, lamu=NULL,
                               lam.min.ratio=1e-8, nlam=10,tol = 1e-9,maxiter=5000){
  return(list(lams=lams, lamu=lamu, 
              lam.min.ratio=lam.min.ratio, nlam=nlam,tol=tol,maxiter=maxiter))
}


#generate lambda sequence
lambda.gen = function(h,q,control.cv=list()){
    lamu1 = svd(matrix(2*h,nrow=prod(q)))$d[1]
    laml1 = lamu1 * control.cv$lam.min.ratio
    alpha = exp(seq(log(1e-04),log(0.5),len=4))
    alpha = c(0,alpha,rev(1-alpha)[-1],1)

    lamseq1 = exp(seq(log(laml1),log(lamu1),len=control.cv$nlam))
    lamseq1 = rev(lamseq1)
    lamseq1[1]=0.8*lamseq1[1]

    lambda=list()
    k=1
    for (i in 1:control.cv$nlam) {
        for (j in 1:length(alpha)) {
            lambda[[k]] = list(lam = lamseq1[i],alpha = alpha[j])
            k=k+1
        }
        alpha = rev(alpha)
    }
    return(lambda)
}




library(foreach)
library(doParallel)
library(Matrix)
library(rTensor)


####################################
##########  5-fold CV  #############
####################################

mOpCov.cv=function(location, x, subject,q,b0 = NULL,ker = "sob", control.cv=list(),control=list()){
  control <- do.call(mOpCov.control, control)
  control.cv <- do.call(mOpCov.cv.control,control.cv)
# Overall 
  prepout0 = prep(location,x,subject,q,Mmethod=control$Mmethod, tol=control$preptol, nsam=control$nsam, ker = ker)
  n = length(prepout0$Ms)
 
# construct folds  
  folds = sample(rep(1:5,length.out=n))
  Testind = list()
  for(i in 1:5){
    Testind[[i]]= which(folds==i)
  }
  if(is.null(b0)) b0 = rep(1e-08,prod(q^2))

lambda = lambda.gen(prepout0$h,q,control.cv)
control$eta = lambda[[1]]$lam*1e-04
#cross validation
registerDoParallel(5)
errormat = foreach(j=1:5,.combine=cbind) %dopar% {
   mOpCov_cv_cpp(prepout0$Xs, prepout0$Ms, q, setdiff(1:n, Testind[[j]]), Testind[[j]], lambda,
    b0, eta=control$eta,
    maxiter=control.cv$maxiter, tol=control.cv$tol)
}
stopImplicitCluster()

error = rowSums(errormat)
lamselect = lambda[[which.min(error)]]
Minv = list()
a = c(0,cumsum(q))
for (i in 1:length(q)) {
  Minv[[i]] =  ginv(prepout0$M[,((a[i]+1):a[i+1])])
}
b0 = rep(1e-08,prod(q^2))

eig = eigen(prepout0$G)
GU =eig$vectors
vG = eig$values
if(control$eta>10*lamselect$lam) {eta = 10*lamselect$lam}else{eta = control$eta}
Final = ACADMM(q,control$maxiter,prepout0$G,GU,vG,prepout0$h,b0,control$eta,lamselect$lam*lamselect$alpha,lamselect$lam*(1-lamselect$alpha)/2,control$tol)

#### further rank reduction for solution
b = Final$b
d = length(q)
perm = c(rev(1:d),rev((1+d):(2*d)))
Bt = as.tensor(array(unfold(b,perm,c(rev(q),rev(q))),c(q,q)))
r = rep(0,d)
for (i in 1:d) {
  Mtemp = k_unfold(Bt,m=i)@data
  r[i] = rankMatrix(Mtemp,tol=1e-5)[1]
}

if(!all(diff(r-q)<0)) {
  for(i in 1:20){
    Bt = as.tensor(array(unfold(b,perm,c(rev(q),rev(q))),c(q,q)))
    Bt_low = tucker(Bt,c(r,r))
    Bm = matrix(unfold(as.vector(Bt_low$est@data),perm,c(q,q)),nrow=prod(q))
    Bm = 0.5 * Bm +0.5 *t(Bm)
    # ensure pd
    eig = eigen(Bm)
    ind = eig$values>=0
    Bm = eig$vectors[,ind] %*% diag(eig$values[ind]) %*% t(eig$vectors[,ind])
    b = as.vector(Bm)
  }
  R = rankMatrix(Bm,tol=1e-10)[1]
  Final$b = b
}

return(list(solution=Final$b,Minv=Minv,q=q,location=prepout0$location,ker = ker, R=R, r =r,lambda_seq=lambda, lambda_select=lamselect,error = error))

}

##################################################
## covariance function prediction on new points ##
##################################################

predict.mOpCov=function(newlocation, OUT){
  d = ncol(newlocation)
  if(OUT$ker == "sob") {ker.fc = K_sob_cpp3}
  if(OUT$ker == "cos") {ker.fc = K_cos_cpp3}
  ktest = sapply(newlocation[,1], ker.fc, OUT$location[,1])
  temp = t(OUT$Minv[[1]] %*% ktest)
  for (i in 2:d) {
    ktest = sapply(newlocation[,i], ker.fc, OUT$location[,i])
    temp = rowkhatrirao(temp,t(OUT$Minv[[i]] %*% ktest))
  }
  Zget = temp %*% matrix(OUT$solution,nrow = prod(OUT$q)) %*% t(temp)
  return(Zget)
}


#########################################################
############ FPCA and marginal basis ####################
#########################################################

## Q matrix 
 K_cos=function(s, t){
  k1 = abs(s-t)/2.0
  k2 = (s+t)/2.0
  s1 = k1^4 - 2* k1^3 + k1^2 - 1.0/30.0
  s2 = k2^4 - 2* k2^3 + k2^2 - 1.0/30.0
  res = (-1.0/3.0)* (s1 + s2)
  return(res)
 }
 
 K_sob=function(s, t){
   k1s = s-0.5
   k1t = t-0.5
   k1abst = abs(s-t)-0.5
   k2s = (k1s * k1s - 1.0/12)/2
   k2t = (k1t * k1t - 1.0/12)/2
   k4abst = (k1abst^4 - k1abst*k1abst/2 + 7.0/240)/24
   return (1.0 + k1s *k1t + k2s *k2t - k4abst)
 }
library(pracma)

 getQ = function(location,  o = 10000, ker, fix.pos=F){
   N = nrow(location)
   tt = seq(0,1,length.out = o)
   Qlist= list()
   for (i in 1:ncol(location)) {
     if(ker == "sob") { Z = outer(location[,i], tt, K_sob)}
     if(ker == "cos") { Z = outer(location[,i], tt, K_cos)}
     #Z = outer(location[,i], tt, K_sob)
     out =  Z %*% t(Z) / o
     if (fix.pos){
       ee <- eigen(out, symmetric=T)
       ee$values[ee$values<0] <- 0
       out <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
     }
     Qlist[[i]]= out
   }
   return(Qlist)
 }

 
 fpca.mOpCov = function(OUT,Qlist=NULL){
   ### eigen functions
   if (is.null(Qlist)) Qlist = getQ(OUT$location,ker = OUT$ker)
   Rkron = OUT$Minv[[1]] %*% Qlist[[1]] %*% t(OUT$Minv[[1]])
   Rlist = list()
   Rlist[[1]] = Rkron
   for (i  in 2:length(OUT$q)) {
     Rlist[[i]] = OUT$Minv[[i]] %*% Qlist[[i]] %*% t(OUT$Minv[[i]])
     Rkron = kronecker(Rkron,Rlist[[i]])
   }
   Rsqr = sqrtm(Rkron)$B
   twom = Rsqr%*% matrix(OUT$solution,nrow = prod(OUT$q)) %*% Rsqr
   e = eigen(twom,symmetric=T)
   rank = OUT$R
   Eigen = list(values = e$values[1:rank], vectors = e$vectors[,1:rank])
   
   #### calculate marginal basis 
   ONE= list()
   d = length(OUT$q)
   perm = c(rev(1:d),rev((1+d):(2*d)))
   for (i in 1:d) {
     ONE_m = matrix(as.vector(aperm(array(as.vector(twom),dim = c(rev(OUT$q),rev(OUT$q))), perm = perm)),OUT$q[i])
     perm_n = perm
     perm_n[1]= perm[i+1]
     perm_n[i+1] = perm[1]
     perm = perm_n
     R1sqr =  sqrtm(Rlist[[i]])$B
     svds = svd(ONE_m)
     ONE[[i]] = list(values = svds$d[1:(OUT$r[i])],vectors = svds$u[,1:(OUT$r[i])])
   }
    return(list(Eigen = Eigen, ONE = ONE, Rlist = Rlist, Rsqr= Rsqr) )
 }
 
 
 #compute eigenfunction (value on newlocation: (Nnew * Nnew))
 computeEigen=function(newlocation,OUT,fpcaOUT){
   if(OUT$ker == "sob") {ker.fc = K_sob_cpp3}
   if(OUT$ker == "cos") {ker.fc = K_cos_cpp3}
   z =  sapply(newlocation[,1], ker.fc, OUT$location[,1]) #(Ndata * Nnew)
   temp = OUT$Minv[[1]] %*% z
   #temp = crossprod(ginv(fpcaOUT$Rsqrlist[[1]]),OUT$Minv[[1]]) %*% z  #(q * Nnew)
   for (i in 2:length(OUT$q)) {
     z = sapply(newlocation[,i], ker.fc, OUT$location[,i])
     temp = khatri_rao(temp,OUT$Minv[[i]] %*% z)
   }
   Res = crossprod(fpcaOUT$Eigen$vectors, ginv(fpcaOUT$Rsqr)%*%temp)
   return(t(Res))
 }
 
#### plot eigenfunctions (two-dimenstional)
 

 
 #### compute marginal basis functions on newlocation
 ## ind is one-way indicator (indicates which dimension)
 computeONE = function(newpoints, OUT, fpcaOUT, ind){
   if(OUT$ker == "sob") {ker.fc = K_sob_cpp3}
   if(OUT$ker == "cos") {ker.fc = K_cos_cpp3}
   u1 = fpcaOUT$ONE[[ind]]$vectors
   z =  sapply(newpoints, ker.fc , OUT$location[,ind])
   temp = OUT$Minv[[ind]] %*% z
   R1sqr =  sqrtm(fpcaOUT$Rlist[[ind]])$B
   ones= crossprod(u1, ginv(R1sqr)%*%temp)
   return(t(ones))
 }
 
 
 
 
