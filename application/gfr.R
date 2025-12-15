#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=gfr
#SBATCH --output=logs/gfr_%j.out
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

library(fda)
library(Matrix)
library(dplyr)
library(tidyr)

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
  choices = c(0,1),
  help = "Whether to run competing methods."
)
args <- parser$parse_args()
compare <- as.logical(args[["compare"]])

cat("Whether compare to OnlineCov:", compare, "\n")



# par.old <- par(no.readonly = TRUE)


## Read GFR data ============
dat <- read.csv("data/gfr/GFR2023Feb01.csv")
# remove the subjects with less than 1 observation
rm_idx <- apply(is.na(dat[, 3:9]), 1, sum) >= 7
dat <- dat[!rm_idx, ]

N <- 143600
# N <- 50000
# N <- 400
m <- 7
t0 <- 0
t1 <- 1
tgrid <- seq(t0, t1, length.out = 7)

# de-mean
gfr_mean <- colMeans(as.matrix(dat[, 3:9]), na.rm = TRUE)
dat$Ly <- unname(apply(
  as.matrix(dat[, 3:9]),
  1,
  \(x) (x - gfr_mean)[!is.na(x)],
  simplify = FALSE
))
dat$Ltid <- unname(apply(
  as.matrix(dat[, 3:9]),
  1,
  \(x) which(!is.na(x)),
  simplify = FALSE
))
dat$Lmi <- sapply(dat$Ltid, length)
# rm(dat)

# barplot(
#   table(dat$Lmi) / length(dat$Lmi) * 100,
#   xlab = "Number of Records",
#   ylab = "Percent of Recipients (%)"
# )
# mean(dat$Lmi == 7)
# mean(dat$Lmi == 6)
# mean(dat$Lmi < 6)

# re-scale
dat$Lt <- lapply(dat$Lt, \(t) (t - 1) / (7 - 1))
ysd <- sd(unlist(dat$Ly))
dat$Ly <- lapply(dat$Ly, \(y) y / ysd)

fdata_generator <- function(n, total_count) {
  idx <- (total_count):(total_count + n - 1) %% N + 1
  return(list(
    Ly = dat$Ly[idx],
    Ltid = dat$Ltid[idx],
    Lt = dat$Lt[idx],
    Lmi = dat$Lmi[idx]
  ))
}

# Algorithm settings =========
nBatch <- 20
nParams <- 3
nPass <- 1
nBlock <- 200
Ninit <- 100
q <- 3
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize0 <- 1e-1
stepsize.min <- 1e-4
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
asgd.use <- TRUE

nBlockIter <- nBlock / nBatch
nIter1pass <- round(N / nBatch)
asgdIterStart <- round(nRecord.1pass * 0.3 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)

# create basis
nbasis <- 7
basis <- create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = 4)
p <- nbasis

B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
GR <- Matrix::chol(G)
invG <- chol2inv(GR)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv

perm <- get_commute_index(p, q)

evalGrid <- seq(t0, t1, length.out = 101)


## Initialization =============
set.seed(123)
idx <- which(dat$Lmi == 7)[1:Ninit]

init_start <- Sys.time()
ppcaObj <- Rdimtools::do.ppca(do.call(rbind, dat$Ly[idx]), ndim = q)
init_end <- Sys.time()
PhiInitEval <- ppcaObj$projection
ThetaInit <- smooth_basis(tgrid, PhiInitEval, basis, lambda = 1e-8)$coefs
lambdaInit <- colMeans(ppcaObj$Y^2)
norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- as.matrix(B %*% ThetaInit)
# matplot(
#   PhiInitEval,
#   type = "l",
#   lty = 3,
#   xlab = "t",
#   ylab = expression(phi[k](t)),
#   main = "Initial Estimates"
# )
sigma2Init <- ppcaObj$mle.sigma2

ThetaInit <- manifold.Stiefel.retract(ThetaInit, NULL, G)
grad_Theta_init <- objfun(
  dat$Ly[1:Ninit], dat$Ltid[1:Ninit],
  ThetaInit, lambdaInit, sigma2Init, NULL, 0, "grad"
)$grad_Theta
grad_Theta_init <- manifold.Stiefel.project(grad_Theta_init, ThetaInit, G)
g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init) / q)
sgd_lr0 <- stepsize0 / g_Theta_init_norm
message("> Step size = ", stepsize0)

## Online FPCA ==============

inits <- list(
  Theta = ThetaInit,
  lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
  sigma2 = sigma2Init
)

message("Start Online FPCA")

fit <- fpca.sgd(
  fdata_generator,
  tgrid,
  inits = inits,
  meanfun = FALSE,
  tau = NULL,
  tau.control = list(ntau = nParams, nselect = 1, maxtau = 1e-0, mintau = 1e-3),
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
  adamw = TRUE,
  adam.rescale = TRUE,
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

resCI <- fit$CI

tau.min <- fit$tau.min
l <- which(fit$tau == tau.min)[1]
Theta.avg <- fit$Theta.avg[,, l]
lambda.avg <- fit$lambda.avg[, l]
sigma2.avg <- fit$sigma2.avg[l]
PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg, basis))

params <- fit$params.history$params
vcrits <- fit$vcrit.history

# check eigenfunctions
tau_path <- with(fit$tau.select, extract_tau_path(tau.history, tau.selectId, l))
tau_path_id <- c(tau_path$tau_path_id, rep(1, (nPass - 1) * round(N / nBlock)))
tau_path_id_extend <- c(rep(tau_path_id[1], nRoundNoTune), tau_path_id)
tau.select = fit$tau.select

ThetaAll.avg <- sapply(
  seq_along(fit$params.history$iter.params),
  \(i) params$Theta.avg[,, tau_path_id_extend[i], i],
  simplify = "array"
)
PhiAvgAll <- array(
  apply(ThetaAll.avg, 3, \(X) eval_fd(evalGrid, FuncData(X, basis))),
  dim = c(length(evalGrid), q, dim(ThetaAll.avg)[3])
)

sgd_time <- colSums(fit$time.history[, 1:nIter1pass])
sgd_time <- colSums(matrix(sgd_time, nrow = nBlockIter))
sgd_time <- c(difftime(init_end, init_start, units = 'secs'), sgd_time)

save(
  PhiAvgEval,
  PhiAvgAll,
  lambda.avg,
  sgd_time,
  tau.select,
  tau_path,
  resCI,
  file = file.path("application", "result_gfr.Rdata")
)


if (compare) {

  ## Online FPCA, local polynomials ------------------------

  message("Start Online Polynomials")

  EV1 <- 100
  EV2 <- 50
  eval_mu <- seq(t0, t1, length.out = EV1) # grid for mean function
  eval_gam_vec <- seq(t0, t1, length.out = EV2)
  eval_gam_mat <- cbind(rep(eval_gam_vec, each = EV2), rep(eval_gam_vec, EV2)) # grid for cov function
  Kmax <- round(N / nBlock) # total number of data blocks, 1000

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
    C = 1.5,
    verbose = TRUE,
    period = 1
  )

  # check eigenfunctions
  PhiEstAll.ll <- sapply(
    1:Kmax,
    \(k) matrix(fit.ll$vecs[1:(q * EV2), k], ncol = q),
    simplify = "array"
  )
  PhiEst.ll <- matrix(fit.ll$vecs[1:(q * EV2), Kmax], ncol = q)
  lambdaEst.ll <- fit.ll$vals[1:q, Kmax]
  PhiEstFunc.ll <- smooth_basis(eval_gam_vec, PhiEst.ll, basis, lambda = 1e-10)
  PhiEst.ll <- eval_fd(evalGrid, PhiEstFunc.ll)

  loclin_time <- fit.ll$time

  # check FPC PVE
  lambdaEst.ll.all = pmax(fit.ll$vals[, Kmax], 0)
  round((lambdaEst.ll.all[1:q] / sum(lambdaEst.ll.all)), 4)
  # 0.912 0.0750 0.0118

  matched_pcs <- match_fpc(PhiEst.ll[, 1:q], PhiAvgEval)
  PhiEst.ll[, 1:q] <- matched_pcs
  PhiEstAll.ll <- PhiEstAll.ll[, attributes(matched_pcs)$match_id, , drop = F]
  lambdaEst.ll[1:q] <- lambdaEst.ll[1:q][attributes(matched_pcs)$match_id]

  save(
    PhiEst.ll,
    PhiEstAll.ll,
    lambdaEst.ll,
    loclin_time,
    file = file.path("application", "result_gfr_olcov.Rdata")
  )

}
