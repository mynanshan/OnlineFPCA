
library(fda)
library(Matrix)


source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./R/onlineFPCA.R")
source("./R/onlineFDAlocalpoly.R")
source("./data_generation/generator.R")

# 1D Simulation -----------------------------------------------------------


dirpath <- file.path("experiments", "example-taupath")

if (!dir.exists(dirpath))
  dir.create(dirpath)

tau_paths = list(
  "1D" = list("RSGD" = NULL, "RAdaGrad" = NULL),
  "2D" = list("RSGD" = NULL, "RAdaGrad" = NULL)
)


## Experiment settings
noise_sd <- 0.5
alpha <- 2
nBatch <- 5
nParams <- 6
nPass <- 3
nBlock <- 100
Ninit <- 100
initMethod <- "face"
q <- 3
N <- 5000
# N <- 2000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize <- 2e-1
stepsize.min <- 5e-2
sgd.step.scale <- 1 # use a smaller step size for sgd
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
asgd.use <- TRUE

message("alpha=", alpha, ", ", "noise=", noise_sd)

set.seed(1)

rp <- get_rp.yang2021(alpha = alpha, npc = 10)
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
fdata_generator <- function(n, total_count) {
  idx <- (total_count):(total_count + n - 1) %% N + 1
  return(list(
    Ly = dat$Ly[idx],
    Ltid = dat$Ltid[idx],
    Lt = dat$Lt[idx],
    Lmi = dat$Lmi[idx]
  ))
}
tgrid <- dat$tgrid

# set B-spline basis
nbasis <- 7
basis <- create.bspline.basis(c(t0[1], t1[1]), nbasis = nbasis, norder = 4)
p <- nbasis

# best possible basis approximations
muTrueFunc <- smooth_basis(evalGrid, muTrueEval, basis)
phiTrueFunc <- smooth_basis(evalGrid, PhiTrueEval, basis)
rmseMuBest <- Metrics::rmse(muTrueEval, eval_fd(evalGrid, muTrueFunc))
rmsePhiBest <- Metrics::rmse(PhiTrueEval, eval_fd(evalGrid, phiTrueFunc))
# Metrics::rmse(PhiTrueEval[,1:q], eval_fd(evalGrid, phiTrueFunc)[,1:q])

# some invariants: basis matrices
B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)

## Initialization ---------------------------------------------------
start_time.init <- Sys.time()
if (initMethod == "face") {
  # Init method 1: FACE
  tmpdat <- data.frame(
    argvals = unlist(dat$Lt[1:Ninit]),
    subj = rep(1:Ninit, dat$Lmi[1:Ninit]),
    y = unlist(dat$Ly[1:Ninit])
  )
  fitFace <- face::face.sparse(tmpdat,
    center = FALSE,
    argvals.new = evalGrid,
    knots = p)
} else if (initMethod == "pace") {
  # Init method 1: PACE
  fitPace <- fdapace::FPCA(dat$Ly[1:Ninit],
    dat$Lt[1:Ninit],
    optns = list(
      dataType = "Sparse",
      nRegGrid = neval,
      userMu = list(t = evalGrid, mu = muTrueEval)
    ))
}
end_time.init <- Sys.time()

## Online FPCA, adaptive smoothness, penalized spline ------------------------
if (initMethod == "face") {
  muInitEval <- fitFace$mu.new
  theta_muInit <- smooth_basis(evalGrid, muInitEval, basis, lambda = 1e-10)$coefs
  eig <- RSpectra::eigs_sym(fitFace$Chat.new, q)
  PhiInitEval <- flip_direc(eig$vectors, PhiTrueEval[, 1:q])
  ThetaInit <- smooth_basis(evalGrid, PhiInitEval, basis, lambda = 1e-10)$coefs
  lambdaInit <- eig$values
  sigma2Init <- fitFace$sigma2
} else if (initMethod == "pace") {
  muInitEval <- fitPace$mu
  theta_muInit <- smooth_basis(evalGrid, muInitEval, basis, lambda = 1e-10)$coefs
  eig <- RSpectra::eigs_sym(fitPace$smoothedCov, q)
  PhiInitEval <- flip_direc(eig$vectors, PhiTrueEval[, 1:q])
  ThetaInit <- smooth_basis(evalGrid, PhiInitEval, basis, lambda = 1e-10)$coefs
  lambdaInit <- eig$values
  sigma2Init <- fitPace$sigma2
}
if (any(lambdaInit <= .Machine$double.eps))
  next

norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- eval_fd(evalGrid, FuncData(ThetaInit, basis))

message("sgd, stepsize=", stepsize, ", ", "nbatch=", nBatch)

nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)

inits <- list(Theta = ThetaInit,
  lambda = lambdaInit,
  sigma2 = sigma2Init)

# Adam --------------
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
  stepsize.decayrate = 0.6,
  stepsize.min = stepsize.min,
  nIter.slowerdecay = floor(0.8 * nIter1pass),
  stepsize.decayrate.slow = 0.3,
  nIter.adam = adamIterEnd,
  asgd.use = TRUE,
  asgd.start = asgdIterStart,
  coord.scaling = FALSE,
  nIter.1stTune = nRoundNoTune * nBlockIter,
  nIter.lastTune = nIter1pass,
  nIter.tauNoIncrease = floor(0.7 * nIter1pass),
  period.tune = nBlockIter,
  period.record = nBlockIter,
  verbose = FALSE
)

tau.min <- fit$tau.min
l <- which(fit$tau == tau.min)[1]

tau_paths$`1D`$RAdam = list(
  tau.history = fit$tau.select$tau.history, 
  tau.selectId = fit$tau.select$tau.selectId,
  final.id = l, nbatch_per_block = nBlock / nBatch
)

setEPS()
postscript(file.path(dirpath, "taupath-sim1d-radam.eps"), width = 5, height = 5)
do.call(plot.tau_path2, tau_paths$`1D`$RAdam)
dev.off()

# SGD ---------------------
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
  stepsize = stepsize * sgd.step.scale,
  stepsize.decayrate = 0.6,
  stepsize.min = stepsize.min,
  nIter.slowerdecay = floor(0.8 * nIter1pass),
  stepsize.decayrate.slow = 0.3,
  nIter.adam = 0,
  asgd.use = TRUE,
  asgd.start = asgdIterStart,
  coord.scaling = FALSE,
  nIter.1stTune = nRoundNoTune * nBlockIter,
  nIter.lastTune = nIter1pass,
  nIter.tauNoIncrease = floor(0.7 * nIter1pass),
  period.tune = nBlockIter,
  period.record = nBlockIter,
  verbose = FALSE
)

tau.min <- fit$tau.min
l <- which(fit$tau == tau.min)[1]

tau_paths$`1D`$RSGD = list(
  tau.history = fit$tau.select$tau.history, 
  tau.selectId = fit$tau.select$tau.selectId,
  final.id = l, nbatch_per_block = nBlock / nBatch
)

do.call(plot.tau_path2, tau_paths$`1D`$RSGD)
