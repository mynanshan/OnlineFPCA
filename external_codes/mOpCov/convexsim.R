Q1=3
Q2=2

library(funData)
library(rTensor)

evl= 1/(1:(Q1*Q2))^2


U1 = diag(Q1)
U2 = diag(Q2)



n=100
m=10
library(MASS)

coreCOE=mvrnorm(n,mu=rep(0,Q1*Q2),Sigma = diag(evl)) #### xi
fullCOE=coreCOE%*%t(kronecker(U1,U2)) #### gamma


S1=NULL
S2=NULL

### rearrange evl to match order of U1 and U2
Ucoe = order(kronecker(1:Q1,1:Q2))
evl = evl[Ucoe]

PHI1=function(v){
    re=NULL
    for (k in 1:length(v)) {
        re=cbind(re,sqrt(2)*cos(pi*v[k]*(1:Q1)))
    }
    return(re)
}

PHI2=function(v){
    re=NULL
    for (k in 1:length(v)) {
        re=cbind(re,sqrt(2)*cos(pi*v[k]*(1:Q2)))
    }
    return(re)
}
#grid

Z=rep(0,n*m)
for (i in 1:n) {
    s1=runif(m)
    s2=runif(m)
    S1=c(S1,s1)
    S2=c(S2,s2)
    coe=khatri_rao(PHI1(s1),PHI2(s2))
    z=colSums(apply(coe,2,function(x) x*fullCOE[i,]))+rnorm(m,0,0.1)
    #z=apply(coe,2,function(x) sum(x*fullCOE[i,]))+rnorm(m,0,0.01) #### alternative
    Z[((i-1)*m+1):(i*m)]=z
}



#for test
s1test=seq(0.02,0.98,length.out = 50)
s2test=seq(0.02,0.98,length.out = 50)
Mph1=PHI1(s1test)
Mph2=PHI2(s2test)
KK=kronecker(t(Mph1),t(Mph2))
Ztest=KK%*%kronecker(U1,U2)%*%diag(evl)%*%t(kronecker(U1,U2))%*%t(KK)


Rcpp::sourceCpp("convex_cpp.cpp")
source("convex_prep.R")



location = cbind(S1,S2)
subject = rep(1:n,each = m)
q = c(6,6)
CVRes = mCovEst.cv(location, Z, subject,q,b0 = NULL)


cat(paste("lam=",CVRes$lambda_select$lam, "\n", sep = ""))
cat(paste("alpha=",CVRes$lambda_select$alpha, "\n", sep = ""))

newlocation = cbind(rep(s1test,each = length(s2test)),rep(s2test,length(s1test)))
Zget = predict.mCovEst(newlocation, CVRes)
DIF=Zget-Ztest
MSE=norm(DIF-diag(diag(DIF)),"F")^2
MSE=MSE/((50^2)*(50^2-1))

cat(paste("MSE_o=",MSE,"\n",sep = ""))

# to get exact low-rank
library(Matrix)
M1=matrix(unfold(CVRes$solution,c(2,1,4,3),c(6,6,6,6)),nrow=q[1])
R1=rankMatrix(M1,tol=1e-10)[1]
M2=matrix(unfold(CVRes$solution,c(1,2,4,3),c(6,6,6,6)),nrow=q[2])
R2=rankMatrix(M2,tol=1e-10)[1]
library(Matrix)
M0=matrix(CVRes$solution,nrow=prod(q))
Rank22=rankMatrix(M0,tol=1e-10)[1]

cat(paste("Rank0_o=",Rank22,"\n",sep = " "))
cat(paste("Rank1_o=",R1,"\n",sep = " "))
cat(paste("Rank2_o=",R2,"\n",sep = " "))


R1=rankMatrix(M1,tol=1e-5)[1]
M2=matrix(unfold(CVRes$solution,c(1,2,4,3),c(6,6,6,6)),nrow=q[2])
R2=rankMatrix(M2,tol=1e-5)[1]

b = CVRes$solution

if(R1<6 | R2<6) {
    for(i in 1:20){
        Bt = as.tensor(array(unfold(b,c(2,1,4,3),c(6,6,6,6)),c(6,6,6,6)))
        Bt_low = tucker(Bt,c(R1,R2,R1,R2))
        Bm = matrix(unfold(as.vector(Bt_low$est@data),c(2,1,4,3),c(6,6,6,6)),nrow=prod(q))
        Bm = 0.5 * Bm +0.5 *t(Bm)
        # ensure pd
        eig = eigen(Bm)
        ind = eig$values>=0
        Bm = eig$vectors[,ind] %*% diag(eig$values[ind]) %*% t(eig$vectors[,ind])
        b = as.vector(Bm)
    }
CVRes$solution = b
}



library(Matrix)
M1=matrix(unfold(CVRes$solution,c(2,1,4,3),c(6,6,6,6)),nrow=q[1])
R1=rankMatrix(M1,tol=1e-10)[1]
M2=matrix(unfold(CVRes$solution,c(1,2,4,3),c(6,6,6,6)),nrow=q[2])
R2=rankMatrix(M2,tol=1e-10)[1]
Zget = predict.mCovEst(newlocation, CVRes)
DIF=Zget-Ztest
MSE=norm(DIF-diag(diag(DIF)),"F")^2
MSE=MSE/((50^2)*(50^2-1))

cat(paste("MSE=",MSE,"\n",sep = ""))
M0=matrix(CVRes$solution,nrow=prod(q))
Rank22=rankMatrix(M0,tol=1e-10)[1]

cat(paste("Rank0=",Rank22,"\n",sep = " "))
cat(paste("Rank1=",R1,"\n",sep = " "))
cat(paste("Rank2=",R2,"\n",sep = " "))


