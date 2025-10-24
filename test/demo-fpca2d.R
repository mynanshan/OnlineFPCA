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

# seed <- 23
seed <- ceiling(runif(1) * 9999)

RhpcBLASctl::blas_set_num_threads(12)

mOpCov_path <- "external_codes/mOpCov/"
source(paste0(mOpCov_path, "mOpCov_prep.R"))
Rcpp::sourceCpp(paste0(mOpCov_path, "mOpCov_cpp.cpp"))

noise_sd <- 0.1
nBatch <- 5
nParams <- 6
nPass <- 2
nBlock <- 100
Ninit <- 100
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize0 <- 0.5
stepsize.min <- 1e-5
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune

exprmt <- "fpca2d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("simu_", exprmt, "_sd", seedtext, ".csv")

res <- data.frame(
  noise = numeric(),
  seed = numeric(),
  Method = character(),
  StepSize = numeric(),
  N = numeric(),
  nBatch = numeric(),
  Time = numeric(),
  RMSEphi1 = numeric(),
  RMSEphi2 = numeric(),
  RMSEphi3 = numeric(),
  RMSEphi1.avg = numeric(),
  RMSEphi2.avg = numeric(),
  RMSEphi3.avg = numeric()
)

set.seed(seed)

rp <- get_rp.wang2020(alpha = 2)
m_min <- NULL
m_max <- NULL
m_mean <- 25
m_sd <- 6
m_type <- "gaussian"

t0 <- rp$t0
t1 <- rp$t1
D <- rp$ndim

# evaluate true mean and eigenfunctions
nevalList <- c(51, 51)
neval <- prod(nevalList)
evalGridList <- lapply(1:rp$ndim, \(d) {
  seq(t0[d], t1[d], length.out = nevalList[d])
})
evalGrid <- margins2grid(evalGridList)
muTrueEval <- rp$meanfun(evalGrid)
PhiTrueEval <- rp$eigfun(evalGrid)
lambdaTrue <- rp$eigval
covTrueEval <- PhiTrueEval %*% diag(lambdaTrue) %*% t(PhiTrueEval)
evalGridListSmall <- lapply(1:rp$ndim, \(d) {
  seq(t0[d], t1[d], length.out = 21)
})
evalGridSmall <- margins2grid(evalGridListSmall)

q <- npc <- 3

# set B-spline basis
nbasis1d <- 7
basis <- TensorBasis(list(
  create.bspline.basis(c(t0[1], t1[1]), nbasis = nbasis1d, norder = 4),
  create.bspline.basis(c(t0[2], t1[2]), nbasis = nbasis1d, norder = 4)
))
p <- attr(basis, "nbasis")

muTrueFunc <- smooth_basis(evalGrid, muTrueEval, basis)
phiTrueFunc <- smooth_basis(evalGrid, PhiTrueEval, basis)
rmseMuBest <- Metrics::rmse(muTrueEval, eval_fd(evalGrid, muTrueFunc))
rmsePhiBest <- Metrics::rmse(PhiTrueEval, eval_fd(evalGrid, phiTrueFunc))

ThetaTrue = phiTrueFunc$coefs[,1:q]
sigma2true = noise_sd^2

G <- get_basis_inprod_matrix(basis)
GR <- chol(G)
invG <- chol2inv(GR)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv

message(">> noise = ", noise_sd)

set.seed(seed)
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
B <- eval_basis(tgrid, basis)

fdata_generator <- function(n, total_count) {
  idx <- (total_count):(total_count + n - 1) %% N + 1
  list(
    Ly = dat$Ly[idx],
    Ltid = dat$Ltid[idx],
    Lt = dat$Lt[idx],
    Lmi = dat$Lmi[idx]
  )
}

# Online FPCA  ------------------------------------------------------------

message("Start OnlineFPCA ---------------")

start_time.init <- Sys.time()
Yinit = dat$Ly[1:Ninit]
Tinit = dat$Lt[1:Ninit]
LmiInit = dat$Lmi[1:Ninit]
TidInit <- dat$Ltid[1:Ninit]
nsub = round((LmiInit + runif(Ninit, min = -0.1, max = 0.1)) * 0.5)
nsub = pmax(nsub, 4)
nsub = pmin(nsub, 20)
nsub = pmin(nsub, LmiInit)
sampleIds = mapply(
  \(mi, mi_sub) {
    sort(sample(1:mi, mi_sub))
  },
  LmiInit,
  nsub
)
Yinit = unlist(mapply(
  \(yi, ids) {
    yi[ids]
  },
  Yinit,
  sampleIds
))
Tinit = do.call(
  rbind,
  mapply(
    \(ti, ids) {
      ti[ids, ]
    },
    Tinit,
    sampleIds
  )
)
TidInit = unlist(mapply(
  \(tid, ids) {
    tid[ids]
  },
  TidInit,
  sampleIds
))
subject_ids <- rep(1:Ninit, times = nsub)
CovRes <- mOpCov(
  location = Tinit,
  x = Yinit,
  subject = subject_ids,
  q = c(6, 6),
  lam = list(lam = 1e-10, alpha = 1e-6),
  ker = "cos"
)
# CovRes <- mOpCov(
#   location = Tinit,
#   x = Yinit,
#   subject = subject_ids,
#   q = c(8, 8),
#   lam = list(lam = 1e-10, alpha = 1e-6),
#   # ker = "cos",
#   ker = "sob",
#   control = list(Mmethod = "eig")
# )
FpcaOut <- fpca.mOpCov(OUT = CovRes)
PhiInitEvalFull <- computeEigen(evalGrid, CovRes, FpcaOut)
lambdaInit <- FpcaOut$Eigen$values[1:q]
PhiInitEval <- PhiInitEvalFull[, 1:q]
ThetaInit <- smooth_basis(
  evalGrid,
  PhiInitEval,
  basis,
  lambda = 1e-8
)$coefs
norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- eval_fd(evalGrid, FuncData(ThetaInit, basis))
PhiInitEval <- flip_direc(PhiInitEval, PhiTrueEval[, 1:q])
flipId <- attributes(PhiInitEval)$flipped
ThetaInit[, flipId] <- -ThetaInit[, flipId]
# fit diagonal element, estimate sigma2
innerIds = which(
  evalGrid[, 1] > 0.1 &
    evalGrid[, 1] < 0.9 &
    evalGrid[, 2] > 0.1 &
    evalGrid[, 2] < 0.9
)
innerPoints = which(TidInit %in% innerIds)
covInitEvalDiag <- as.vector(PhiInitEvalFull^2 %*% FpcaOut$Eigen$values)
sigma2Init <- max(mean((Yinit^2 - covInitEvalDiag[TidInit])[innerPoints]), 1e-3)
end_time.init <- Sys.time()

norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- eval_fd(evalGrid, FuncData(ThetaInit, basis))

nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch

ThetaInit <- manifold.Stiefel.retract(ThetaInit, NULL, G)
grad_Theta_init <- objfun(
  dat$Ly[1:Ninit], dat$Ltid[1:Ninit],
  ThetaInit, lambdaInit, sigma2Init, NULL, 0, "grad"
)$grad_Theta
grad_Theta_init <- manifold.Stiefel.project(grad_Theta_init, ThetaInit, G)
g_Theta_init_norm <- sqrt(mean(colSums(grad_Theta_init * G %*% grad_Theta_init)))
# g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init))
sgd_lr0 <- stepsize0 / g_Theta_init_norm
message("> Step size = ", stepsize0)

# TODO: it turns out the g_Theta_init_norm is not a good estimate
# of the grad's norm
# Also, the algorithm is not necessaily better when ||direc|| is approximatedly 1
# TODO: check an adaptive base step size scheme

inits <- list(
  Theta = ThetaInit,
  lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
  # sigma2 = min(sigma2Init, 1e-2)
  sigma2 = sigma2Init
)
# inits <- list(
#   Theta = ThetaTrue[,1:q] + rnorm(p*q, 0, 1e-2),
#   lambda = lambdaTrue[1:q] * exp(rnorm(q, 0, 1e-2)),
#   sigma2 = sigma2true * exp(rnorm(1, 0, 1e-2))
# )

sgdtype <- "sgd"
sgdtype <- "adam"
sgdtype <- "adagrad"

message(">>> sgdtype = ", sgdtype)

fit <- fpca.sgd(
  fdata_generator,
  tgrid,
  inits = inits,
  meanfun = FALSE,
  tau = NULL,
  tau.control = list(
    ntau = nParams,
    nselect = 2,
    maxtau = 1e-3,
    mintau = 1e-8
  ),
  nbatch = nBatch,
  maxIter = nPass * nIter1pass,
  stepsize = if(sgdtype %in% c("sgd", "sgdm")) {sgd_lr0} else { # stepsize0
    0.1
  },
  nIter.constStepSize = 0,
  stepsize.decayrate = 0.51, # 0.51,
  stepsize.min = stepsize.min,
  period.decay = 5 * nBlockIter,
  nIter.slowerdecay = nIter1pass, # 10 * nBlockIter, # same as adareset
  stepsize.decayrate.slow = 0.3,
  dynlr = ifelse(sgdtype %in% c("sgd", "sgdm"), TRUE, FALSE),
  # dynlr = FALSE,
  dynlrCtrl = list(
    niter = 20 * nBlockIter,
    reset = 5 * nBlockIter,
    refdn = 0.1,
    w = 0.9
  ),
  sgdtype = sgdtype,
  adamw = TRUE,
  adam.rescale = TRUE,
  adareset = 10 * nBlockIter,
  adareset.end = nIter1pass,
  asgd.start = 10 * nBlockIter,
  asgd.reset = 20 * nBlockIter,
  asgd.reset.end = nIter1pass,
  asgd.end = Inf,
  # weight = "subj",
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
rmseAll <- rmse_phi(ThetaAll, phiTrueFunc$coefs, B)
ThetaAll.avg <- sapply(
  seq_along(fit$params.history$iter.params),
  \(i) params$Theta.avg[,, tau_path_id_extend[i], i],
  simplify = "array"
)
rmseAll.avg <- rmse_phi(ThetaAll.avg, phiTrueFunc$coefs, B)

matplot(rmseAll, type="l")
matplot(rmseAll.avg, type="l")

par(mfrow=c(1,q))
for (k in seq_len(q)) {
  plot(PhiAvgEval[,k], xlab="", ylab=bquote(phi[.(k)]), col = rgb(0, 0, 0, 0.4))
  points(PhiTrueEval[,k], col = rgb(1, 0, 0, 0.4))
}
par(mfrow=c(1,1))

# (ThetaAll[,,1:150] - ThetaAll[,,2:151]) |> 
#   matrix(nrow = p) |> 
#   (\(C) as.matrix(B %*% C)^2)() |> 
#   colMeans() |> sqrt() |> 
#   matrix(nrow = q) |> t()

for (ii in seq(51, 501, 25)) {
k <- 2
plot(as.vector(B %*% ThetaAll[,k,ii]), xlab="", ylab=bquote(phi[.(k)]), col = rgb(0, 0, 0, 0.4))
points(PhiTrueEval[,k], col = rgb(1, 0, 0, 0.4))
}



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
  dat$Ly, dat$Ltid, phiTrueFunc$coefs, lambdaTrue, noise_sd^2,
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
    Theta[,2] <- Theta[,2] + h * Delta1
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

objfun(
  dat$Ly, dat$Ltid, ThetaBatch0, lambdaTrue[1:q], noise_sd^2,
  tau = tau.min, stats = "loss"
)$lik