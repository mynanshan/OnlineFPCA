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

noise_sd <- 0.5
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
stepsize0 <- 0.1
stepsize.min <- 1e-4
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune

exprmt <- "ci1d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("ci_", exprmt, "_sd", seedtext, ".RData")

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

inits <- list(
  Theta = ThetaInit,
  lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
  sigma2 = sigma2Init
)

sgd_lr0 <- stepsize0 / g_Theta_init_norm
message("> Step size = ", stepsize0)

sgdtype <- "sgd"
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
  # stepsize = ifelse(sgdtype=="sgd", sgd_lr0, stepsize0),
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
  asgd.reset = Inf,
  asgd.reset.end = Inf,
  asgd.end = Inf,
  fpcCI = TRUE,
  nIter.1stTune = nRoundNoTune * nBlockIter,
  nIter.lastTune = nIter1pass,
  nIter.tauNoIncrease = floor(0.3 * nIter1pass),
  dyntune.rate = 1.25,
  period.tune = nBlockIter,
  period.record = nBlockIter,
  verbose = TRUE
)

resCI <- fit$CI

tau.min <- fit$tau.min
l <- which(fit$tau == tau.min)[1]
Theta.avg <- fit$Theta.avg[,, l]

# match eigenfunctions
PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg, basis))
PhiAvgEval <- match_fpc(PhiAvgEval, PhiTrueEval[, 1:q, drop = F])
resCI <- resCI[,
  attributes(PhiAvgEval)$match_id, , , drop = F ]
resCI[, attributes(PhiAvgEval)$flipped, , ] <-
  -resCI[, attributes(PhiAvgEval)$flipped, , ]


# ---------- Compute a global optimum ----------------
library(ManifoldOptim)

Th_id <- seq(1, p*q)
lam_id <- seq(p*q + 1, p*q + q)
sig_id <- p*q + q + 1

F_val <- function(x) {
  Theta <- as.matrix(backsolve(GR, matrix(x[Th_id], p, q)))
  lambda <- exp(x[lam_id])
  sigma2 <- exp(x[sig_id])
  objfun(dat$Ly, dat$Ltid, Theta, lambda, sigma2, NULL, 1e-9, "loss")$fval
}
F_grad <- function(x) {
  Theta <- backsolve(GR, matrix(x[Th_id], p, q))
  lambda <- exp(x[lam_id])
  sigma2 <- exp(x[sig_id])
  gradObj <- objfun(dat$Ly, dat$Ltid, Theta, lambda, sigma2, NULL, 1e-9, "grad")
  c(
    as.matrix(GR %*% gradObj$grad_Theta),
    gradObj$grad_eta,
    gradObj$grad_zeta
  )
}
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, F_val, F_grad)

mani.defn <- get.product.defn(
  get.stiefel.defn(p, q), get.euclidean.defn(q+1, 1)
)
mani.params <- get.manifold.params()
solver.params <- get.solver.params(isconvex = FALSE, DEBUG=3)
x0 <- c(
  as.matrix(GR %*% phiTrueFunc$coefs[,1:q]),
  log(lambdaTrue[1:q]), log(noise_sd^2)
)
opt <- manifold.optim(
  prob, mani.defn, method = "LRBFGS", 
  mani.params = mani.params, solver.params = solver.params, x0 = x0
)

ThetaStar <- as.matrix(backsolve(GR, matrix(opt$xopt[Th_id], p, q)))
lambdaStar <- exp(opt$xopt[lam_id])
sigma2Star <- exp(opt$xopt[sig_id])

# matplot(evalGrid, PhiTrueEval[,1:q], type="l", lty=1)
# matplot(evalGrid, eval_fd(evalGrid, FuncData(ThetaStar, basis)), type="l", add=T, lty=2)

# End experiment ----------------------------------------------------------

CIobj = list(
  PhiStar = as.matrix(B %*% ThetaStar),
  PhiTrue = as.matrix(B %*% phiTrueFunc$coefs[,1:npc]),
  CIs = resCI
) 

save(CIobj, file = file.path(dirpath, filename))

message("Finishing replication: ", seed)
set.seed(NULL)
