#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=news
#SBATCH --output=logs/news_%j.out
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

suppressPackageStartupMessages({
  library(fda)
  library(Matrix)
})

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./R/onlineFPCA.R")
source("./R/onlineFDAlocalpoly.R")

parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)
parser$add_argument(
  "--compare",
  type = "integer",
  default = 1,
  choices = c(0, 1),
  help = "Whether to run competing methods."
)
args <- parser$parse_args()
compare <- as.logical(args[["compare"]])

cat("Whether compare to OnlineCov:", compare, "\n")

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ================================================================
## User settings
## ================================================================
project_root <- "."
processed_file <- file.path(project_root, "data", "news", "processed", "news_preprocessed.rds")
out_dir <- file.path(project_root, "application")
use_partial_last_batch <- TRUE

## OnlineFPCA settings
nBatch <- 5
nBlock <- 50           # reporting/tuning block size (subjects)
nParams <- 3
nPass <- 1
Ninit <- 200
q <- 4                 # Yang–Yao show 4 leading eigenfunctions
nbasis <- 11
basis_order <- 4
stepsize0 <- 1e-1
stepsize.min <- 1e-4
nRoundNoTune <- 1

## Initialization grid used only for PPCA warm start
n_init_grid <- 49

## ================================================================
## Helpers
## ================================================================
subset_stream_dat <- function(dat, idx) {
  list(
    Ly = dat$Ly[idx],
    Ltid = dat$Ltid[idx],
    Lt = dat$Lt[idx],
    Lmi = dat$Lmi[idx]
  )
}

make_generator <- function(dat, Nmax) {
  force(dat)
  force(Nmax)
  function(n, total_count) {
    idx <- seq.int(total_count, min(total_count + n - 1, Nmax))
    list(
      Ly = dat$Ly[idx],
      Ltid = dat$Ltid[idx],
      Lt = dat$Lt[idx],
      Lmi = dat$Lmi[idx]
    )
  }
}

make_init_matrix <- function(dat, idx, xout) {
  t(vapply(idx, function(i) {
    approx(
      x = dat$Lt[[i]],
      y = dat$Ly[[i]],
      xout = xout,
      method = "linear",
      rule = 2
    )$y
  }, numeric(length(xout))))
}

aggregate_times <- function(time_vec, block_iter) {
  gid <- ceiling(seq_along(time_vec) / block_iter)
  as.numeric(tapply(time_vec, gid, sum))
}

match_history_length <- function(tau_path_id, n_slices, nRoundNoTune) {
  if (length(tau_path_id) == 0) {
    tau_path_id <- 1
  }
  out <- c(rep(tau_path_id[1], nRoundNoTune), tau_path_id)
  if (length(out) < n_slices) {
    out <- c(out, rep(tail(out, 1), n_slices - length(out)))
  }
  out[seq_len(n_slices)]
}

## ================================================================
## Load preprocessed data
## ================================================================
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
news_obj <- readRDS(processed_file)
dat_all <- news_obj$dat
N_total <- length(dat_all$Ly)

if (!use_partial_last_batch) {
  N <- floor(N_total / nBatch) * nBatch
} else {
  N <- N_total
}

if (N < q + 5) {
  stop("Too few subjects after preprocessing.")
}

if (N < N_total) {
  message("Dropping the last ", N_total - N, " subjects to keep full mini-batches.")
}

dat <- subset_stream_dat(dat_all, seq_len(N))
fdata_generator <- make_generator(dat, N)

message("Loaded ", N_total, " preprocessed subjects.")
message("Using ", N, " subjects in the online run.")
message("Block size for SGD mini-batches: ", nBatch)
message("Reporting/tuning block size: ", nBlock)

## ================================================================
## Algorithm settings
## ================================================================
t0 <- news_obj$t0
t1 <- news_obj$t1
tgrid <- news_obj$tgrid
nt <- length(tgrid)

nBlockIter <- nBlock / nBatch
if (abs(nBlockIter - round(nBlockIter)) > 1e-8) {
  stop("`nBlock` must be a multiple of `nBatch`.")
}
nBlockIter <- as.integer(round(nBlockIter))

nIter1pass <- ceiling(N / nBatch)
nIters <- seq_len(nIter1pass)

basis <- create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = basis_order)
p <- nbasis
B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
GR <- Matrix::chol(G)
invG <- chol2inv(GR)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv

evalGrid <- seq(t0, t1, length.out = 241)
clock_eval <- 6 + 24 * evalGrid

## ================================================================
## Initialization
## ================================================================
set.seed(123)
init_candidates <- which(dat$Lmi >= 20)
if (length(init_candidates) < Ninit) {
  Ninit <- length(init_candidates)
}
init_idx <- init_candidates[seq_len(Ninit)]
tobs_init <- seq(t0, t1, length.out = n_init_grid)

init_start <- Sys.time()
Y_init <- make_init_matrix(dat, init_idx, xout = tobs_init)
ppcaObj <- Rdimtools::do.ppca(Y_init, ndim = q)
init_end <- Sys.time()

PhiInitEval <- ppcaObj$projection
ThetaInit <- smooth_basis(tobs_init, PhiInitEval, basis, lambda = 1e-8)$coefs
lambdaInit <- colMeans(ppcaObj$Y^2)
norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
sigma2Init <- max(ppcaObj$mle.sigma2, 1e-6)

ThetaInit <- manifold.Stiefel.retract(ThetaInit, NULL, G)
grad_Theta_init <- objfun(
  dat$Ly[init_idx], dat$Ltid[init_idx],
  ThetaInit, lambdaInit, sigma2Init, NULL, 0, "grad"
)$grad_Theta
grad_Theta_init <- manifold.Stiefel.project(grad_Theta_init, ThetaInit, G)
g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init) / q)
sgd_lr0 <- stepsize0 / max(g_Theta_init_norm, 1e-8)

inits <- list(
  Theta = ThetaInit,
  lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
  sigma2 = sigma2Init
)

message("Initialization done. Effective initial learning rate = ", signif(sgd_lr0, 4))

## ================================================================
## Online FPCA without CI
## ================================================================
message("Start OnlineFPCA-RSGD on news data")
fit <- fpca.sgd(
  fdata_generator,
  tgrid,
  inits = inits,
  meanfun = FALSE,
  tau = NULL,
  tau.control = list(ntau = nParams, nselect = 1, maxtau = 1e+0, mintau = 1e-3),
  nbatch = nBatch,
  maxIter = nPass * nIter1pass,
  stepsize = sgd_lr0,
  nIter.constStepSize = 0,
  stepsize.decayrate = 0.51,
  stepsize.min = stepsize.min,
  period.decay = 5 * nBlockIter,
  nIter.slowerdecay = nIter1pass,
  stepsize.decayrate.slow = 0.25,
  dynlr = TRUE,
  dynlrCtrl = list(
    niter = 20 * nBlockIter,
    reset = 5 * nBlockIter,
    refdn = stepsize0,
    w = 0.9
  ),
  sgdtype = "sgd",
  ada.start = 25 * nBlockIter + 1,
  adareset = 20 * nBlockIter,
  adareset.end = nIter1pass,
  asgd.start = 10 * nBlockIter + 1,
  asgd.reset = Inf,
  asgd.reset.end = nIter1pass,
  asgd.end = Inf,
  fpcCI = FALSE,
  nIter.1stTune = nRoundNoTune * nBlockIter,
  nIter.lastTune = nIter1pass,
  nIter.tauNoIncrease = floor(1 * nIter1pass),
  period.tune = nBlockIter,
  period.record = nBlockIter,
  verbose = TRUE
)

tau.min <- fit$tau.min
l <- which(fit$tau == tau.min)[1]
Theta.avg <- fit$Theta.avg[, , l]
lambda.avg <- fit$lambda.avg[, l]
sigma2.avg <- fit$sigma2.avg[l]
PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg, basis))

params <- fit$params.history$params
tau.select <- fit$tau.select
tau_path <- with(tau.select, extract_tau_path(tau.history, tau.selectId, l))

n_slices <- length(fit$params.history$iter.params)
tau_path_id_extend <- match_history_length(tau_path$tau_path_id, n_slices, nRoundNoTune)
ThetaAll.avg <- sapply(
  seq_len(n_slices),
  function(i) params$Theta.avg[, , tau_path_id_extend[i], i],
  simplify = "array"
)
PhiAvgAll <- array(
  apply(ThetaAll.avg, 3, function(X) eval_fd(evalGrid, FuncData(X, basis))),
  dim = c(length(evalGrid), q, dim(ThetaAll.avg)[3])
)

iter_time <- colSums(fit$time.history[, seq_len(nIter1pass), drop = FALSE])
sgd_time <- c(as.numeric(difftime(init_end, init_start, units = "secs")), aggregate_times(iter_time, nBlockIter))

record_subjects <- if (dim(PhiAvgAll)[3] == length(sgd_time)) {
  c(0, pmin(seq_len(length(sgd_time) - 1L) * nBlock, N))
} else if (dim(PhiAvgAll)[3] == length(sgd_time) - 1L) {
  pmin(seq_len(dim(PhiAvgAll)[3]) * nBlock, N)
} else {
  pmin(seq_len(dim(PhiAvgAll)[3]) * nBlock, N)
}

run_config <- list(
  processed_file = processed_file,
  N_total = N_total,
  N = N,
  nBatch = nBatch,
  nBlock = nBlock,
  nPass = nPass,
  q = q,
  nbasis = nbasis,
  basis_order = basis_order,
  nBlockIter = nBlockIter,
  nIter1pass = nIter1pass,
  evalGrid = evalGrid,
  clock_eval = clock_eval,
  tobs_init = tobs_init,
  use_partial_last_batch = use_partial_last_batch
)

save(
  PhiAvgEval,
  PhiAvgAll,
  lambda.avg,
  sigma2.avg,
  sgd_time,
  tau.min,
  tau_path,
  tau.select,
  record_subjects,
  run_config,
  file = file.path(out_dir, "result_news.Rdata")
)

## ================================================================
## Online FPCA with CI
## ================================================================
message("Start OnlineFPCA-RSGD with pointwise CI")
fit_uq <- fpca.sgd(
  fdata_generator,
  tgrid,
  inits = inits,
  meanfun = FALSE,
  tau = NULL,
  tau.control = list(ntau = nParams, nselect = 1, maxtau = 1e+0, mintau = 1e-3),
  nbatch = nBatch,
  maxIter = nPass * nIter1pass,
  stepsize = sgd_lr0,
  nIter.constStepSize = 0,
  stepsize.decayrate = 0.51,
  stepsize.min = stepsize.min,
  period.decay = 5 * nBlockIter,
  nIter.slowerdecay = nIter1pass,
  stepsize.decayrate.slow = 0.25,
  dynlr = TRUE,
  dynlrCtrl = list(
    niter = 20 * nBlockIter,
    reset = 5 * nBlockIter,
    refdn = stepsize0,
    w = 0.9
  ),
  sgdtype = "sgd",
  ada.start = 25 * nBlockIter + 1,
  adareset = 20 * nBlockIter,
  adareset.end = nIter1pass,
  asgd.start = 10 * nBlockIter + 1,
  asgd.reset = Inf,
  asgd.reset.end = nIter1pass,
  asgd.end = Inf,
  fpcCI = TRUE,
  nIter.1stTune = nRoundNoTune * nBlockIter,
  nIter.lastTune = nIter1pass,
  nIter.tauNoIncrease = floor(1 * nIter1pass),
  dyntune.rate = 1.25,
  period.tune = nBlockIter,
  period.record = nBlockIter,
  verbose = FALSE
)

resCI <- fit_uq$CI

tau.min.uq <- fit_uq$tau.min
l.uq <- which(fit_uq$tau == tau.min.uq)[1]
Theta.avg.uq <- fit_uq$Theta.avg[, , l.uq]
lambda.avg.uq <- fit_uq$lambda.avg[, l.uq]
PhiAvgEval.uq <- eval_fd(evalGrid, FuncData(Theta.avg.uq, basis))

params.uq <- fit_uq$params.history$params
tau_path.uq <- with(fit_uq$tau.select, extract_tau_path(tau.history, tau.selectId, l.uq))
n_slices_uq <- length(fit_uq$params.history$iter.params)
tau_path_id_extend.uq <- match_history_length(tau_path.uq$tau_path_id, n_slices_uq, nRoundNoTune)
ThetaAll.avg.uq <- sapply(
  seq_len(n_slices_uq),
  function(i) params.uq$Theta.avg[, , tau_path_id_extend.uq[i], i],
  simplify = "array"
)
PhiAvgAll.uq <- array(
  apply(ThetaAll.avg.uq, 3, function(X) eval_fd(evalGrid, FuncData(X, basis))),
  dim = c(length(evalGrid), q, dim(ThetaAll.avg.uq)[3])
)

iter_time.uq <- colSums(fit_uq$time.history[, seq_len(nIter1pass), drop = FALSE])
sgd_time.uq <- c(as.numeric(difftime(init_end, init_start, units = "secs")), aggregate_times(iter_time.uq, nBlockIter))
record_subjects.uq <- if (dim(PhiAvgAll.uq)[3] == length(sgd_time.uq)) {
  c(0, pmin(seq_len(length(sgd_time.uq) - 1L) * nBlock, N))
} else if (dim(PhiAvgAll.uq)[3] == length(sgd_time.uq) - 1L) {
  pmin(seq_len(dim(PhiAvgAll.uq)[3]) * nBlock, N)
} else {
  pmin(seq_len(dim(PhiAvgAll.uq)[3]) * nBlock, N)
}

tau.select.uq <- fit_uq$tau.select

save(
  PhiAvgEval.uq,
  PhiAvgAll.uq,
  lambda.avg.uq,
  sgd_time.uq,
  tau.min.uq,
  tau_path.uq,
  tau.select.uq,
  resCI,
  record_subjects.uq,
  run_config,
  file = file.path(out_dir, "result_news_uq.Rdata")
)

## ================================================================
## Optional comparison: OnlineCov / local polynomial covariance
## ================================================================
if (compare) {
  message("Start OnlineCov comparison on centered/scaled data")

  EV1 <- 121
  EV2 <- 61
  eval_mu <- seq(t0, t1, length.out = EV1)
  eval_gam_vec <- seq(t0, t1, length.out = EV2)
  eval_gam_mat <- cbind(rep(eval_gam_vec, each = EV2), rep(eval_gam_vec, EV2))
  Kmax <- ceiling(N / nBlock)

  fit.ll <- fpca.lpoly.online(
    fdata_generator,
    n = nBlock,
    evalArgs = list(
      EV1 = EV1,
      EV2 = EV2,
      t0 = t0,
      t1 = t1,
      eval_mu = eval_mu,
      eval_gam_vec = eval_gam_vec,
      eval_gam_mat = eval_gam_mat
    ),
    streamArgs = list(Kmax = Kmax),
    L2 = nParams,
    G = 0.5,
    R = 1,
    Mcl = 1,
    C0 = 1.5,
    verbose = TRUE,
    period = 1
  )

  PhiEstAll.ll <- sapply(
    seq_len(Kmax),
    function(k) matrix(fit.ll$vecs[1:(q * EV2), k], ncol = q),
    simplify = "array"
  )
  PhiEst.ll <- matrix(fit.ll$vecs[1:(q * EV2), Kmax], ncol = q)
  lambdaEst.ll <- fit.ll$vals[1:q, Kmax]
  PhiEstFunc.ll <- smooth_basis(eval_gam_vec, PhiEst.ll, basis, lambda = 1e-10)
  PhiEst.ll <- eval_fd(evalGrid, PhiEstFunc.ll)
  loclin_time <- fit.ll$time

  matched_pcs <- match_fpc(PhiEst.ll[, 1:q], PhiAvgEval)
  PhiEst.ll[, 1:q] <- matched_pcs
  PhiEstAll.ll <- PhiEstAll.ll[, attributes(matched_pcs)$match_id, , drop = FALSE]
  lambdaEst.ll[1:q] <- lambdaEst.ll[1:q][attributes(matched_pcs)$match_id]

  record_subjects.ll <- pmin(seq_len(Kmax) * nBlock, N)

  save(
    PhiEst.ll,
    PhiEstAll.ll,
    lambdaEst.ll,
    loclin_time,
    record_subjects.ll,
    run_config,
    file = file.path(out_dir, "result_news_olcov.Rdata")
  )
}

message("Done.")
