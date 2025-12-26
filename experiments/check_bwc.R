#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=test-bwc
#SBATCH --output=logs/test-bwc_%j.out
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G

parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)
parser$add_argument(
  "--seed",
  type = "integer",
  default = 1234,
  help = "Seed used for this experiment."
)
args <- parser$parse_args()
seed <- as.integer(args[["seed"]])

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

nBatch <- 5
noise_sd <- 0.5
alpha <- 2
nParams <- 6
nPass <- 3
nBlock <- 100
Ninit <- 100
initMethod <- "face"
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize0 <- 1e-1
stepsize.min <- 1e-4
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune

bwcSet <- data.frame(
  C = c(3, 4, 4, 6, 6, 6, 8, 8, 8),
  W = c(1, 1, 2, 1, 2, 3, 1, 2, 4)
)

exprmt <- "bwc"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

if (seed == 1) {
  readr::write_csv(bwcSet, file.path(dirpath, "bwcSet.csv"))
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
basename <- paste0(exprmt, "_sd", seedtext)

res <- data.frame(
  seed = numeric(),
  Method = character(),
  StepSize = numeric(),
  N = numeric(),
  Ninit = numeric(),
  initMethod = character(),
  npc = numeric(),
  nBatch = numeric(),
  Time = numeric(),
  C = integer(),
  W = integer(),
  RMSEphi1 = numeric(),
  RMSEphi2 = numeric(),
  RMSEphi3 = numeric(),
  RMSEphi1.avg = numeric(),
  RMSEphi2.avg = numeric(),
  RMSEphi3.avg = numeric()
)

tau_paths_mat <- c()

set.seed(seed)

rp <- get_rp.yang2021(alpha = alpha, npc = 10)
m_min <- NULL
m_max <- NULL
m_mean <- 6
m_sd <- 2
m_type <- "gaussian"

t0 <- rp$t0
t1 <- rp$t1
D <- rp$ndim

q <- 3

# evaluate true mean and eigenfunctions
neval <- 101
evalGrid <- seq(t0, t1, length.out = neval)
muTrueEval <- rp$meanfun(evalGrid)
PhiTrueEval <- rp$eigfun(evalGrid)
lambdaTrue <- rp$eigval

# set B-spline basis
nbasis <- 7
basis <- create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = 4)
p <- nbasis

muTrueFunc <- smooth_basis(evalGrid, muTrueEval, basis)
phiTrueFunc <- smooth_basis(evalGrid, PhiTrueEval, basis)
rmseMuBest <- Metrics::rmse(muTrueEval, eval_fd(evalGrid, muTrueFunc))
rmsePhiBest <- Metrics::rmse(PhiTrueEval, eval_fd(evalGrid, phiTrueFunc))

G <- get_basis_inprod_matrix(basis)
GR <- Matrix::chol(G)
invG <- chol2inv(GR)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv

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
B <- eval_basis(tgrid, basis)


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
theta_muInit <- smooth_basis(evalGrid, muInitEval, basis, lambda = 1e-10)$coefs
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
g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init) / q)
sgd_lr0 <- stepsize0 / g_Theta_init_norm
message("> Step size = ", stepsize0)

nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch

inits <- list(
  Theta = ThetaInit,
  lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
  sigma2 = sigma2Init
)

for (sgdtype in c("sgd", "adagrad")) {
  message("Method:", sgdtype)

  for (i in seq_len(nrow(bwcSet))) {
    message(">>>", "C:", bwcSet$C[i], ", W:", bwcSet$W[i])

    fit <- fpca.sgd(
      fdata_generator,
      tgrid,
      inits = inits,
      meanfun = FALSE,
      tau = NULL,
      tau.control = list(
        ntau = bwcSet$C[i],
        nselect = bwcSet$W[i],
        maxtau = 1e-1,
        mintau = 1e-6
      ),
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
      sgdtype = sgdtype,
      adamw = TRUE,
      adam.rescale = TRUE,
      ada.start = 25 * nBlockIter + 1,
      adareset = 20 * nBlockIter,
      adareset.end = nIter1pass,
      asgd.start = 10 * nBlockIter + 1,
      asgd.reset = 40 * nBlockIter,
      asgd.reset.end = nIter1pass,
      asgd.end = Inf,
      nIter.1stTune = nRoundNoTune * nBlockIter,
      nIter.lastTune = nIter1pass,
      nIter.tauNoIncrease = floor(0.3 * nIter1pass),
      period.tune = nBlockIter,
      period.record = nBlockIter,
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
    params$Theta <- params$Theta[, attributes(PhiEstEval)$match_id, , , drop = F]
    params$Theta[, attributes(PhiEstEval)$flipped, , ] <-
      -params$Theta[, attributes(PhiEstEval)$flipped, , ]
    PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg, basis))
    PhiAvgEval <- match_fpc(PhiAvgEval, PhiTrueEval[, 1:q, drop = F])
    params$Theta.avg <- params$Theta.avg[, attributes(PhiAvgEval)$match_id, , , drop = F]
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

    times <- c(
      colSums(fit$time.history[, 1:nIter1pass]),
      fit$time.history[1, (nIter1pass + 1):(nPass * nIter1pass)]
    )
    times <- colSums(matrix(times, nrow = nBlockIter))
    times <- c(difftime(end_time.init, start_time.init, units = 'secs'), times)

    res <- rbind(
      res,
      data.frame(
        seed = rep(seed, nRecord + 1),
        Method = rep(paste0("OnlineFPCA-", sgdtype), nRecord + 1),
        StepSize = rep(stepsize0, nRecord + 1),
        N = c(0, nIters),
        Ninit = rep(Ninit, nRecord + 1),
        initMethod = rep(initMethod, nRecord + 1),
        npc = rep(q, nRecord + 1),
        nBatch = rep(nBatch, nRecord + 1),
        Time = times,
        C = bwcSet$C[i],
        W = bwcSet$W[i],
        RMSEphi1 = rmseAll[, 1],
        RMSEphi2 = rmseAll[, 2],
        RMSEphi3 = rmseAll[, 3],
        RMSEphi1.avg = rmseAll.avg[, 1],
        RMSEphi2.avg = rmseAll.avg[, 2],
        RMSEphi3.avg = rmseAll.avg[, 3]
      )
    )

    tau_paths_mat <- rbind(tau_paths_mat, tau_path)
    rownames(tau_paths_mat)[nrow(tau_paths_mat)] <- paste0(sgdtype, "-C", bwcSet$C[i], "W", bwcSet$W[i])
  }
}

# End experiment ----------------------------------------------------------

readr::write_csv(res, file.path(dirpath, paste0("simu_", basename, ".csv")))
saveRDS(tau_paths_mat, file.path(dirpath, paste0("taupath_", basename, ".csv")))

message("Finishing replication: ", seed)
set.seed(NULL)
