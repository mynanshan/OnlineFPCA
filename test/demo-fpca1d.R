library(fda)
library(Matrix)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./R/onlineFPCA.R")
source("./R/onlineFDAlocalpoly.R")
source("./R/fpcaReg.R")
source("./data_generation/generator.R")

# TODO:
# run the new AdaGrad algo
# Are there closed-form solutions for lambda and sigma2?

noise_sd <- 0.5
nBatch <- 5
nParams <- 6
nPass <- 5
nBlock <- 100
Ninit <- 100
initMethod <- "face"
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize <- 1e-1
stepsize.min <- 1e-2
sgd.step.scale <- 1 # use a smaller step size for sgd
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
asgd.use <- TRUE

# seed <- 2501
# set.seed(seed)

rp <- get_rp.yang2021(alpha = 2, npc = 10)
m_min <- NULL
m_max <- NULL
m_mean <- 6
m_sd <- 2
m_type <- "gaussian"

t0 <- rp$t0
t1 <- rp$t1
D <- rp$ndim

# evaluate true mean and eigenfunctions
neval <- 101
evalGrid <- seq(t0, t1, length.out = neval)
muTrueEval <- rp$meanfun(evalGrid)
PhiTrueEval <- rp$eigfun(evalGrid)
lambdaTrue <- rp$eigval

q <- npc <- 3

dat <- get_measurements(
  rp,
  n = N,
  m_min = m_min,
  m_max = m_max,
  m_mean = m_mean,
  m_sd = m_sd,
  design_type = "random",
  m_type = m_type,
  sigma = noise_sd
)
tgrid <- dat$tgrid

fdata_generator <- function(n, total_count) {
  idx <- (total_count):(total_count + n - 1) %% N + 1
  list(
    Ly = dat$Ly[idx],
    Ltid = dat$Ltid[idx],
    Lt = dat$Lt[idx],
    Lmi = dat$Lmi[idx]
  )
}

nbasis <- 7
basis <- create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = 4)
p <- nbasis

muTrueFunc <- smooth_basis(evalGrid, muTrueEval, basis)
phiTrueFunc <- smooth_basis(evalGrid, PhiTrueEval, basis)
rmseMuBest <- Metrics::rmse(muTrueEval, eval_fd(evalGrid, muTrueFunc))
rmsePhiBest <- sqrt(colMeans(PhiTrueEval[,1:q] - eval_fd(evalGrid, phiTrueFunc)[,1:q])^2)

B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv


tmpdat <- data.frame(
  argvals = unlist(dat$Lt[1:Ninit]),
  subj = rep(1:Ninit, dat$Lmi[1:Ninit]),
  y = unlist(dat$Ly[1:Ninit])
)
fitFace <- face::face.sparse(
  tmpdat,
  center = FALSE,
  argvals.new = evalGrid,
  knots = p
)

muInitEval <- fitFace$mu.new
theta_muInit <- smooth_basis(
  evalGrid,
  muInitEval,
  basis,
  lambda = 1e-10
)$coefs
eig <- RSpectra::eigs_sym(fitFace$Chat.new, q)
PhiInitEval <- flip_direc(eig$vectors, PhiTrueEval[, 1:q])
ThetaInit <- smooth_basis(evalGrid, PhiInitEval, basis, lambda = 1e-10)$coefs
lambdaInit <- eig$values
sigma2Init <- fitFace$sigma2

norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- eval_fd(evalGrid, FuncData(ThetaInit, basis))

nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)

inits <- list(Theta = ThetaInit, lambda = lambdaInit, sigma2 = sigma2Init)

optMethod <- "SGD"

fit <- fpca.sgd(
  fdata_generator,
  tgrid,
  inits = inits,
  meanfun = FALSE,
  tau = NULL,
  tau.control = list(
    ntau = nParams,
    nselect = 2,
    maxtau = 1e-1,
    mintau = 1e-6
  ),
  nbatch = nBatch,
  maxIter = nPass * nIter1pass,
  stepsize = stepsize,
  stepsize.decayrate = 0.5,
  stepsize.min = stepsize.min,
  nIter.slowerdecay = floor(0.8 * nIter1pass),
  stepsize.decayrate.slow = 0.3,
  nIter.adam = ifelse(optMethod == "SGD", 0, adamIterEnd),
  asgd.use = TRUE,
  asgd.start = asgdIterStart,
  coord.scaling = TRUE,
  nIter.1stTune = nRoundNoTune * nBlockIter,
  nIter.lastTune = nIter1pass,
  nIter.tauNoIncrease = floor(0.7 * nIter1pass),
  period.tune = nBlockIter,
  period.record = nBlockIter,
  verbose = TRUE
)

tau.min <- fit$tau.min
l <- which(fit$tau == tau.min)[1]
Theta <- fit$Theta[,, l]
lambda <- fit$lambda[, l]
sigma2 <- fit$sigma2[l]
Theta.avg <- fit$Theta.avg[,, l]
lambda.avg <- fit$lambda.avg[, l]
sigma2.avg <- fit$sigma2.avg[l]

params <- fit$params.history$params
vcrits <- fit$vcrit.history

# check eigenfunctions
PhiEstEval <- eval_fd(evalGrid, FuncData(Theta, basis))
PhiEstEval <- match_fpc(PhiEstEval, PhiTrueEval[, 1:q, drop = F])
params$Theta <- params$Theta[,
  attributes(PhiEstEval)$match_id, , , drop = F ]
params$Theta[, attributes(PhiEstEval)$flipped, , ] <-
  -params$Theta[, attributes(PhiEstEval)$flipped, , ]
PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg, basis))
PhiAvgEval <- match_fpc(PhiAvgEval, PhiTrueEval[, 1:q, drop = F])
params$Theta.avg <- params$Theta.avg[,
  attributes(PhiAvgEval)$match_id, , , drop = F ]
params$Theta.avg[, attributes(PhiAvgEval)$flipped, , ] <-
  -params$Theta.avg[, attributes(PhiAvgEval)$flipped, , ]

tau_path <- with(
  fit$tau.select,
  extract_tau_path(tau.history, tau.selectId, l)
)
tau_path_id <- c(
  tau_path$tau_path_id,
  rep(1, (nPass - 1) * round(N / nBlock))
)
tau_path_id_extend <- c(rep(tau_path_id[1], nRoundNoTune), tau_path_id)

ThetaAll <- sapply(
  seq_along(fit$params.history$iter.params),
  \(i) params$Theta[,, tau_path_id_extend[i], i],
  simplify = "array"
)
rmseAll <- sapply(1:q, \(k) {
  diffTheta <- sweep(ThetaAll[, k, ], 1, phiTrueFunc$coefs[, k], "-")
  diffPhi <- B %*% diffTheta
  sqrt(colMeans(diffPhi^2))
})
ThetaAll.avg <- sapply(
  seq_along(fit$params.history$iter.params),
  \(i) params$Theta.avg[,, tau_path_id_extend[i], i],
  simplify = "array"
)
rmseAll.avg <- sapply(1:q, \(k) {
  diffTheta <- sweep(ThetaAll.avg[, k, ], 1, phiTrueFunc$coefs[, k], "-")
  diffPhi <- B %*% diffTheta
  sqrt(colMeans(diffPhi^2))
})

matplot(evalGrid, PhiTrueEval[,1:q], type="l", lty=1)
matplot(tgrid, eval_fd(tgrid, FuncData(Theta.avg, basis)), type="l", add=T, lty=2)
matplot(rmseAll, type="l")
matplot(rmseAll.avg, type="l")

# check opt landscape
objfun(
  dat$Ly, dat$Ltid, Theta.avg, lambda.avg, sigma2.avg,
  tau = tau.min, stats = "loss"
)
objfun(
  dat$Ly, dat$Ltid, ThetaBatch, lambda.avg, sigma2.avg,
  tau = tau.min, stats = "loss"
)
objfun(
  dat$Ly, dat$Ltid, Theta.avg, lambda.avg, sigma2.avg,
  tau = tau.min, stats = "grad"
)

tmp <- mapply(
  \(h) objfun(
    dat$Ly, dat$Ltid, Theta.avg, lambda.avg, sigma2.avg * exp(h),
    tau = tau.min, stats = "loss"
  )$fval,
  h = seq(-0.2,0.2,0.02)
)



library(rstiefel)

sol <- optStiefel(
  F = \(X) {
    ThetaTmp <- sqrtGinv %*% X
    objfun(
      dat$Ly, dat$Ltid, ThetaTmp, lambda.avg, sigma2.avg,
      tau = tau.min, stats = "loss"
    )$fval
  },
  dF = \(X) {
    ThetaTmp <- sqrtGinv %*% X
    gradObj <- objfun(
      dat$Ly, dat$Ltid, ThetaTmp, lambda.avg, sigma2.avg,
      tau = tau.min, stats = "grad"
    )
    as.matrix(sqrtG %*% gradObj$grad_Theta)
  },
  Vinit = sqrtG %*% ThetaInit,
  verbose = TRUE
)

ThetaBatch <- sqrtGinv %*% sol
matplot(evalGrid, PhiTrueEval[,1:q], type="l", lty=1)
matplot(tgrid, eval_fd(tgrid, FuncData(ThetaBatch, basis)), type="l", add=T, lty=2)
sqrt(colMeans((PhiTrueEval[,1:q] - eval_fd(evalGrid, FuncData(ThetaBatch, basis)))^2))


sol0 <- optStiefel(
  F = \(X) {
    ThetaTmp <- sqrtGinv %*% X
    objfun(
      dat$Ly, dat$Ltid, ThetaTmp, lambdaTrue[1:q], noise_sd^2,
      tau = tau.min, stats = "loss"
    )$fval
  },
  dF = \(X) {
    ThetaTmp <- sqrtGinv %*% X
    gradObj <- objfun(
      dat$Ly, dat$Ltid, ThetaTmp, lambdaTrue[1:q], noise_sd^2,
      tau = tau.min, stats = "grad"
    )
    as.matrix(sqrtG %*% gradObj$grad_Theta)
  },
  Vinit = sqrtG %*% ThetaInit,
  verbose = TRUE
)

ThetaBatch0 <- sqrtGinv %*% sol0
matplot(evalGrid, PhiTrueEval[,1:q], type="l", lty=1)
matplot(tgrid, eval_fd(tgrid, FuncData(ThetaBatch0, basis)), type="l", add=T, lty=2)
sqrt(colMeans((PhiTrueEval[,1:q] - eval_fd(evalGrid, FuncData(ThetaBatch0, basis)))^2))
