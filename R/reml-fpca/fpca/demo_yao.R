setwd("D:/ResearchesAndProjects/2021_9_FPCA/ref_codes/Peng-reml-fpca/fpca")


source("R/fpca_func.R")
source("R/fpca_pred.R")
source("R/fpca_score.R")
source("R/functions_EM.R")
source("R/functions_GenData.R")
source("R/functions_LocLin.R")
source("R/functions_Optimization.R")


load("data/yao2005rescale.Rdata")
dat <- as.matrix(dat)


##sample trajectory
plot(dat[,3],dat[,2],xlab="grid",ylab="",type='p')

for(i in 1:50){
  cur <- dat[dat[,1]==i,]
  points(cur[,3],cur[,2],type="l")
}

## candidate models for fitting
M.set<-c(6,7)
r.set<-c(2)

##parameters for fpca.mle
ini.method="loc"
basis.method="bs"
sl.v=rep(0.5,3)
max.step=50
grid.l=seq(0,1,0.01)
grids=seq(0,1,0.002)


##fit candidate models by fpca.mle
result<-fpca.mle(dat, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
summary(result)

##rescaled grid
grids.new <- result$grid

##model selection result
M<-result$selected_model[1]
r<-result$selected_model[2]


##compare estimated eigenvalues with the truth
(evalest<-result$eigenvalues)     ## estimated
gen$eigvals    ## true

##compare estimated error variance with the truth
(sig2est<-result$error_var)        ## estimated
gen$sigma^2     ## true

##plot: compare estimated eigenfunctions with the truth
eigenfest<-result$eigenfunctions
eigenf<-t(do.call(cbind, lapply(gen$phi, \(f) f(grids.new))))  ##true
for (j in 1:r) {
  if (sum(eigenf[j,] * eigenfest[j,]) < 0) {
    eigenfest[j,] <- -eigenfest[j,]
  }
}
par(mfrow=c(1,2))
for(i in 1:r){
  plot(grids.new,eigenfest[i,],type="l",
       ylim=range(eigenfest),xlab="time",ylab=paste("eigenfunction",i),
       lwd = 2)
  lines(grids, eigenf[i,],col=2,lwd=3,lty=3)
}
par(mfrow=c(1,1))

##plot: compare estimated mean curve with the truth

muest<-result$fitted_mean
plot(grids.new,muest,type="l",lwd=2,
     xlab="time",ylab="mean curve",ylim=range(result$fitted_mean))
lines(grids,numeric(length(grids)),col=5)
par(mfrow=c(1,1))

##look at the CV scores and convergence for each model: note that model (M=5, r=4) does not converge. 
result$cv_scores   ##CV
result$converge   ##convergence


## derive fpc scores and look at the predicted curve
#fpc scores
fpcs<-fpca.score(dat,grids.new,muest,evalest,eigenfest,sig2est,r)
#get predicted trajectories on a fine grid: the same grid for which mean and eigenfunctions are evaluated
pred<-fpca.pred(fpcs, muest,eigenfest)

#get predicted trajectories on the observed measurement points
N<-length(grids.new)

dev.new()
par(mfrow=c(3,3))
for (i in 1:9){
  id<-i                                      ##for curve i
  t.c<-easy$data[easy$data[,1]==id,3]    ##measurement points
  t.proj<-ceiling(N*t.c)                     ##measurement points projected on the grid
  y.c<-easy$data[easy$data[,1]==id,2]    ##obs
  y.pred.proj<-pred[t.proj,id]               ##predicted obs on the measurement points
  
  #plots
  plot(t.c,y.c,ylim=range(pred[,id]),xlab="time",ylab="obs", main=paste("predicted trajectory of curve", id))
  points(grids.new,pred[,id],col=3,type='l')
  ##points(t.c,y.pred.proj,col=2, pch=2)     ##predicted measurements at observed measurement times
}
par(mfrow=c(1,1))
dev.off()

