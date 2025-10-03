#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=check_abv_score_1d
#SBATCH --output=logs/abv1d_%j.out
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G

parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)
parser$add_argument(
  "--seed",
  type = "integer",
  default = 1,
  help = "Seed used for this experiment."
)
args <- parser$parse_args()
seed <- as.integer(args[["seed"]])

library(fda)
library(Matrix)
library(foreach)
library(doFuture)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./R/onlineFPCA.R")
source("./R/onlineFDAlocalpoly.R")
source("./R/fpcaReg.R")
source("./data_generation/generator.R")

noise_sd <- 0.5
nBatch <- 5
nParams <- 6
nPass <- 3
nBlockList <- c(10, 20, 50, 100)
Ninit <- 100
initMethod <- "face"
N <- 5000
stepsize0 <- 0.5
stepsize.min <- 1e-3
ewmabv.beta.list <- c(0.1, 0.3, 0.5, 0.7, 0.9)

exprmt <- "abv1d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("simu_", exprmt, "_sd", seedtext, ".csv")

set.seed(seed)

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
rmsePhiBest <- Metrics::rmse(PhiTrueEval, eval_fd(evalGrid, phiTrueFunc))

B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
GR <- Matrix::chol(G)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv

# Online FPCA  ------------------------------------------------------------

message("Start OnlineFPCA ---------------")

start_time.init <- Sys.time()
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
end_time.init <- Sys.time()

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

ThetaInit <- manifold.Stiefel.retract(ThetaInit, NULL, G)
grad_Theta_init <- objfun(
  dat$Ly[1:Ninit], dat$Ltid[1:Ninit],
  ThetaInit, lambdaInit, sigma2Init, NULL, 0, "grad"
)$grad_Theta
grad_Theta_init <- manifold.Stiefel.project(grad_Theta_init, ThetaInit, G)
g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init))
sgd_lr0 <- stepsize0 / g_Theta_init_norm
message("> Step size = ", stepsize0)

inits <- list(
  Theta = ThetaInit,
  lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
  sigma2 = sigma2Init
)

plan(multisession, workers = availableCores() - 1)

res <-
  foreach(nBlock = nBlockList, .combine = "rbind") %:%
  foreach(ewmabv.beta = ewmabv.beta.list, .combine = "rbind") %dofuture% {

    cat("Settings: ", "nBlock =", nBlock, " ewmabv.beta =", ewmabv.beta, "\n")

    nIters.1pass <- seq(nBlock, N, nBlock)
    nIters <- seq(nBlock, nPass * N, nBlock)
    nRecord.1pass <- round(N / nBlock)
    nRecord <- length(nIters)
    nRoundNoTune <- 1
    nRoundTune <- nRecord.1pass - nRoundNoTune
    nBlockIter <- nBlock / nBatch
    nIter1pass <- N / nBatch

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
      stepsize = stepsize0,
      stepsize.decayrate = 0.51,
      stepsize.min = stepsize.min,
      nIter.slowerdecay = floor(0.8 * nIter1pass),
      stepsize.decayrate.slow = 0.51,
      sgdtype = "adagrad",
      adareset = 20 * nBlockIter,
      adareset.end = Inf,
      asgd.start = 10 * nBlockIter,
      asgd.reset = 20 * nBlockIter,
      asgd.reset.end = nIter1pass,
      asgd.end = Inf,
      nIter.1stTune = nRoundNoTune * nBlockIter,
      nIter.lastTune = nIter1pass,
      nIter.tauNoIncrease = floor(0.3 * nIter1pass),
      period.tune = nBlockIter,
      period.record = nBlockIter,
      ewmabv.beta = ewmabv.beta,
      verbose = FALSE
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
    rmseAll <- rmse_phi(ThetaAll, phiTrueFunc$coefs, B)
    ThetaAll.avg <- sapply(
      seq_along(fit$params.history$iter.params),
      \(i) params$Theta.avg[,, tau_path_id_extend[i], i],
      simplify = "array"
    )
    rmseAll.avg <- rmse_phi(ThetaAll.avg, phiTrueFunc$coefs, B)

    times <- c(
      colSums(fit$time.history[, 1:nIter1pass]),
      fit$time.history[1, (nIter1pass + 1):(nPass * nIter1pass)]
    )
    times <- colSums(matrix(times, nrow = nBlockIter))
    times <- c(difftime(end_time.init, start_time.init, units = 'secs'), times)

    Nlist = c(0, nIters)
    idx = Nlist %in% (N * seq_len(nPass))

    data.frame(
      seed = seed,
      Method = paste0("OnlineFPCA-adagrad"),
      epoch = seq_len(nPass),
      nBlock = nBlock,
      wBeta = ewmabv.beta,
      tau = tau.min,
      RMSEphi1 = rmseAll[idx, 1],
      RMSEphi2 = rmseAll[idx, 2],
      RMSEphi3 = rmseAll[idx, 3],
      RMSEphi1.avg = rmseAll.avg[idx, 1],
      RMSEphi2.avg = rmseAll.avg[idx, 2],
      RMSEphi3.avg = rmseAll.avg[idx, 3]
    )
  }

if (!interactive()) {
  readr::write_csv(res, file.path(dirpath, filename))
}