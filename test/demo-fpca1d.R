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


set.seed(1)

noise_sd <- 0.1
nBatch <- 5
nParams <- 6
nPass <- 2
nBlock <- 100
Ninit <- 100
initMethod <- "face"
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize0 <- 0.1
stepsize.min <- 1e-5
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
GR <- Matrix::chol(G)
invG <- chol2inv(GR)
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

# TODO: update these codes to formal versions
# - Retract ThetaInit
# - Project grad_Theta_init
# - RMSE evaluation, pmin version
# - Remove all unused arguments
# - Set ada and asgd related parameters
# - Change lambdaInit: rather large and inaccurate than small
ThetaInit <- manifold.Stiefel.retract(ThetaInit, NULL, G)
grad_Theta_init <- objfun(
  dat$Ly[1:Ninit], dat$Ltid[1:Ninit],
  ThetaInit, lambdaInit, sigma2Init, NULL, 0, "grad"
)$grad_Theta
grad_Theta_init <- manifold.Stiefel.project(grad_Theta_init, ThetaInit, G)
g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init) / q)
sgd_lr0 <- stepsize0 / g_Theta_init_norm

nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)

# inits <- list(Theta = ThetaInit, lambda = lambdaInit, sigma2 = sigma2Init)
inits <- list(
  Theta = ThetaInit,
  lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
  sigma2 = sigma2Init
)

# TODO: how to account for the covariance between FPCs and eigvals
# TODO: bad initialization of lambda significantly affect the performance
#   does adagrad alleviate this issue?

sgdtype <- "sgd"
sgdtype <- "adam"

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
  stepsize = ifelse(sgdtype %in% c("sgd", "sgdm"), sgd_lr0, stepsize0),
  # stepsize = sgd_lr0,
  nIter.constStepSize = 0,
  stepsize.decayrate = 0.51,
  stepsize.min = stepsize.min,
  period.decay = 5 * nBlockIter,
  nIter.slowerdecay = nIter1pass,
  stepsize.decayrate.slow = 0.25,
  dynlr = ifelse(sgdtype %in% c("sgd", "sgdm"), TRUE, FALSE),
  # dynlr = TRUE,
  dynlrCtrl = list(
    niter = 20 * nBlockIter,
    reset = 5 * nBlockIter,
    refdn = stepsize0,
    w = 0.9
  ),
  sgdtype = sgdtype,
  adamw = TRUE,
  adam.rescale = TRUE,
  # ada.start = 5 * nBlockIter,
  ada.start = 20 * nBlockIter + 1,
  adareset = 10 * nBlockIter,
  # adareset.end = nIter1pass,
  adareset.end = 10 * nBlockIter,
  asgd.start = 10 * nBlockIter + 1,
  asgd.reset = 20 * nBlockIter,
  asgd.reset.end = nIter1pass,
  asgd.end = Inf,
  nIter.1stTune = nRoundNoTune * nBlockIter,
  nIter.lastTune = nIter1pass,
  nIter.tauNoIncrease = floor(0.3 * nIter1pass),
  period.tune = nBlockIter,
  period.record = nBlockIter,
  verbose = TRUE
)

check <- fit$check

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
# ThetaAll <- params$Theta[,, 1,]
rmseAll <- rmse_phi(ThetaAll, phiTrueFunc$coefs, B)
ThetaAll.avg <- sapply(
  seq_along(fit$params.history$iter.params),
  \(i) params$Theta.avg[,, tau_path_id_extend[i], i],
  simplify = "array"
)
# ThetaAll.avg <- params$Theta.avg[,, 1,]
rmseAll.avg <- rmse_phi(ThetaAll.avg, phiTrueFunc$coefs, B)

matplot(evalGrid, PhiTrueEval[,1:q], type="l", lty=1)
matplot(evalGrid, PhiAvgEval, type="l", add=T, lty=2)
matplot(rmseAll, type="l")
matplot(rmseAll.avg, type="l")

# check opt landscape
ThetaTrue = phiTrueFunc$coefs[,1:q]
lambdaTrue
sigma2true = noise_sd^2

objfun(
  dat$Ly, dat$Ltid, Theta.avg, lambda.avg, sigma2.avg,
  tau = tau.min, stats = "loss"
)$lik
objfun(
  dat$Ly, dat$Ltid, ThetaTrue, lambda.avg, sigma2.avg,
  tau = tau.min, stats = "loss"
)$lik
objfun(
  dat$Ly, dat$Ltid, ThetaTrue, lambdaTrue[1:q], noise_sd^2,
  tau = tau.min, stats = "loss"
)$lik
objfun(
  dat$Ly, dat$Ltid, ThetaBatch0, lambdaTrue[1:q], noise_sd^2,
  tau = tau.min, stats = "loss"
)$lik

tmp <- mapply(
  \(h) objfun(
    dat$Ly, dat$Ltid, Theta.avg, lambda.avg, sigma2.avg * exp(h),
    tau = tau.min, stats = "loss"
  )$fval,
  h = seq(-0.2,0.2,0.02)
)
plot(tmp)

Delta1 <- rnorm(p)
tmp <- mapply(
  \(h) {
    Theta <- Theta.avg
    Theta[,1] <- Theta[,1] + h * Delta1
    Theta <- manifold.Stiefel.retract(Theta, NULL, G)
    objfun(
    dat$Ly, dat$Ltid, Theta, lambda.avg, sigma2.avg,
    tau = tau.min, stats = "loss"
  )$fval
  },
  h = seq(-0.2,0.2,0.02)
)
plot(tmp)



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
