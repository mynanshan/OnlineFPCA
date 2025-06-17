setwd("D:/ResearchesAndProjects/2021_9_FPCA/ref_codes/Wang_Wong-LowRankCov-MultiDimFunData")

source("example.R")
source("mOpCov_prep.R")
Rcpp::sourceCpp("mOpCov_cpp.cpp")

#####################################################
#### generate a two dimensional functional data ######
#####################################################

###### same as simulation seting 1 #########


# n: number of samples
# m: number of observations per sample
# sigma: standard deviation of noise

set.seed(35373)
Data = generate.data(n = 100, m = 10, sigma = 0.1)


#############################################
############ Implement mOpCov ###############
#############################################

# location : observed location matrix ((n*m) times dimension)
# x : a vector for observed value on each location ((n*m))
# subject: a vector indicating which sample do observations belong to ((n*m))
# q: a vector setting the number of original computing basis functions 
# lam: a list for tunning parameter lambda and alpha, alpha is within [0,1]
# ker: kernel type

CovRes = mOpCov(location=Data$location, x=Data$Z, subject=Data$subject, q=c(6,6), lam = list(lam = 2.8e-06, alpha = 1e-06), ker = "cos")

## Look at the rank estimations
# two-way rank
CovRes$R
# one-way ranks
CovRes$r



###########################################################
## evaluate the fitted covariance function on new points ##
###########################################################
new_points = cbind(runif(10), runif(10))

# obtain a covariance matrix repect to these new points  
CovEst = predict.mOpCov(newlocation = new_points, OUT = CovRes)


############################################################
#### L2 Eigen functions from fitted covariance function ####
############################################################
FpcaOut = fpca.mOpCov(OUT = CovRes)

## eigen values
eigenvalues = FpcaOut$Eigen$values

## evaluate eigenfunctions
new_points = cbind(rep(seq(0,1,length.out = 50),50),rep(seq(0,1,length.out = 50),each=50))
FpcaEval = computeEigen(new_points,CovRes,FpcaOut)

## plot the first and second eigenfunctions
library(ggplot2)
Eigens = data.frame(dim1 = new_points[,1], dim2 = new_points[,2], values = FpcaEval[,1:2])
plot_Eig1 = ggplot() +
geom_raster(data = Eigens, aes(x = dim1, y = dim2, fill = values.1)) +
scale_fill_gradient(low = "yellow", high = "red") 
plot_Eig1

plot_Eig2 = ggplot() +
  geom_raster(data = Eigens, aes(x = dim1, y = dim2, fill = values.2)) +
  scale_fill_gradient(low = "yellow", high = "red") 
plot_Eig2


##############################################
######## marginal basis function #############
##############################################

## singular values for marginal basis function
# dimension 1
FpcaOut$ONE[[1]]$values
# dimension 2
FpcaOut$ONE[[2]]$values


## evaluation
new_points = seq(0,1,length.out = 100)
# dimension 1
OneEval_1 = computeONE(new_points, CovRes, FpcaOut, 1)
# dimension 2
OneEval_2 = computeONE(new_points, CovRes, FpcaOut, 2)

## plot L2 marginal basis functions
# dimension 1
one_1 = data.frame(x = rep(new_points,dim(OneEval_1)[2]), y = as.vector(OneEval_1), group = as.factor(rep(1:dim(OneEval_1)[2],each = 100)))
plot_one_1 <- ggplot(one_1, aes(x=x, y = y, color = group)) + geom_line()
plot_one_1
# dimension 2
one_2 = data.frame(x = rep(new_points,dim(OneEval_2)[2]), y = as.vector(OneEval_2), group = as.factor(rep(1:dim(OneEval_1)[2],each = 100)))
plot_one_2 <- ggplot(one_2, aes(x=x, y = y, color = group)) + geom_line()
plot_one_2
