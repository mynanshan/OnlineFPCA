#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=ci-simu1d
#SBATCH --output=logs/ci-simu1d_%j.out
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

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

noiseList <- c(0.1, 0.5, 1.)
nBatch <- 5
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
stepsize0 <- 0.5
stepsize.min <- 1e-3
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune

exprmt <- "ci1d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filebasename <- paste0("ci_", exprmt, "_sd", seedtext)

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

nbasis <- 7
basis <- fda::create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = 4)
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

perm <- get_commute_index(p, q)


for (noise_sd in noiseList) {

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

  # Online FPCA  ------------------------------------------------------

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

  nBlockIter <- nBlock / nBatch
  nIter1pass <- N / nBatch

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

  for (sgdtype in c("sgd", "adagrad", "adam")) {
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
        maxtau = 1e-1,
        mintau = 1e-6
      ),
      nbatch = nBatch,
      maxIter = nPass * nIter1pass,
      stepsize = ifelse(sgdtype=="sgd", sgd_lr0, stepsize0),
      stepsize.decayrate = 0.51,
      stepsize.min = stepsize.min,
      nIter.slowerdecay = nIter1pass,
      stepsize.decayrate.slow = 0.51,
      sgdtype = sgdtype,
      adamw = TRUE,
      adam.rescale = TRUE,
      adareset = 20 * nBlockIter,
      adareset.end = nIter1pass,
      asgd.start = 10 * nBlockIter,
      asgd.reset = 20 * nBlockIter,
      asgd.reset.end = nIter1pass,
      asgd.end = Inf,
      fpcCI = TRUE,
      ci.start = 10 * nBlockIter,
      nIter.1stTune = nRoundNoTune * nBlockIter,
      nIter.lastTune = nIter1pass,
      nIter.tauNoIncrease = floor(0.3 * nIter1pass),
      period.tune = nBlockIter,
      period.record = nBlockIter,
      verbose = TRUE
    )

    CIs <- fit$CI

    saveRDS(CIs, file.path(dirpath, paste0(filebasename, "_", sgdtype, ".rds")))

  }

}


# End experiment ----------------------------------------------------------


message("Finishing replication: ", seed)
set.seed(NULL)
