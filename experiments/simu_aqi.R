#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=simu_aqi        # Job name
#SBATCH --output=logs/simu_aqi_%j.out    # Standard output file (%j expands to jobID)
#SBATCH --time=3:00:00                   # Maximum runtime (hh:mm:ss)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G                 # Memory allocation

parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)
parser$add_argument(
  "--seed",
  type = "integer",
  default = 5,
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

# load mOpCov codes
mOpCov_path <- "external_codes/mOpCov/"
source(paste0(mOpCov_path, "mOpCov_prep.R"))
Rcpp::sourceCpp(paste0(mOpCov_path, "mOpCov_cpp.cpp"))

# settings
noise <- 0.1
nBatch <- 5
nParams <- 6
nPass <- 3
nBlock <- 100
Ninit <- 500
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize <- 1e-1
stepsize.min <- 1e-3
sgd.step.scale <- 1 # when < 1, use a smaller step size for sgd
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
asgd.use <- TRUE

exprmt <- "aqi"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("simu_", exprmt, "_sd", seedtext, ".csv")

res <- data.frame(
  seed = numeric(),
  Method = character(),
  StepSize = numeric(),
  N = numeric(),
  nBatch = numeric(),
  Time = numeric(),
  # AV = numeric(),
  # EWMABV = numeric(),
  RMSEphi1 = numeric(),
  RMSEphi2 = numeric(),
  RMSEphi3 = numeric(),
  RMSEphi1.avg = numeric(),
  RMSEphi2.avg = numeric(),
  RMSEphi3.avg = numeric()
)

set.seed(seed)

dat <- gendata.aqi(N)
tgrid <- evalGrid <- dat$tgrid
t0 <- dat$t0
t1 <- dat$t1
neval <- nrow(evalGrid)
muTrueEval <- numeric(neval)
PhiTrueEval <- dat$Phi
lambdaTrue <- dat$lam
q <- npc <- ncol(PhiTrueEval)

fdata_generator <- function(n, total_count) {
  idx <- (total_count):(total_count + n - 1) %% N + 1
  list(
    Ly = dat$Ly[idx],
    Ltid = dat$Ltid[idx],
    Lt = dat$Lt[idx],
    Lmi = dat$Lmi[idx]
  )
}

basis <- TensorBasis(list(
  create.bspline.basis(c(0, 1), nbasis = 6, norder = 4),
  create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
))
p <- attr(basis, "nbasis")

ThetaRecord <- rep(list(list(SGD = NULL, Adam = NULL)), nPass)

# best possible basis approximations
muTrueFunc <- smooth_basis(evalGrid, muTrueEval, basis)
phiTrueFunc <- smooth_basis(evalGrid, PhiTrueEval, basis)
rmseMuBest <- Metrics::rmse(muTrueEval, eval_fd(evalGrid, muTrueFunc))
rmsePhiBest <- Metrics::rmse(PhiTrueEval, eval_fd(evalGrid, phiTrueFunc))

# some invariants: basis matrices
B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
RG <- Matrix::chol(G) # t(RG) %*% RG == G
invRG <- Matrix::solve(RG)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)


# Online FPCA  ------------------------------------------------------------
message("Start OnlineFPCA ---------------")

## Initialization ---------------------------------------------------
message("mOpCov Initialization ...")
start_time.init <- Sys.time()
Yinit = dat$Ly[1:Ninit]
Tinit = dat$Lt[1:Ninit]
LmiInit = dat$Lmi[1:Ninit]
TidInit <- dat$Ltid[1:Ninit]
# do a subsampling within each subject to
# accelerate mOpCov
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
# Yinit <- unlist(dat$Ly[1:Ninit])
# Tinit <- do.call(rbind, dat$Lt[1:Ninit])
# LtIdInit <- unlist(dat$Ltid[1:Ninit])
# subject_ids <- rep(1:Ninit, times = sapply(dat$Ly[1:Ninit], length))
CovRes <- mOpCov(
  location = Tinit,
  x = Yinit,
  subject = subject_ids,
  q = c(8, 8),
  lam = list(lam = 1e-10, alpha = 1e-6),
  # ker = "cos",
  ker = "sob",
  control = list(Mmethod = "eig")
)
FpcaOut <- fpca.mOpCov(OUT = CovRes)
lambdaInit <- FpcaOut$Eigen$values[1:q]
PhiInitEvalFull <- computeEigen(evalGrid, CovRes, FpcaOut)
PhiInitEval <- PhiInitEvalFull[, 1:q]
ThetaInit <- smooth_basis(evalGrid, PhiInitEval, basis, lambda = 1e-8)$coefs
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
sigma2Init <- mean((Yinit^2 - covInitEvalDiag[TidInit])[innerPoints])
sigma2Init <- max(sigma2Init, 1e-3)
end_time.init <- Sys.time()

inits <- list(Theta = ThetaInit, lambda = lambdaInit, sigma2 = sigma2Init)
nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)

for (optMethod in c("SGD", "Adam")) {
  message("Method:", optMethod)

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
      mintau = 1e-8
    ),
    nbatch = nBatch,
    maxIter = nPass * nIter1pass,
    stepsize = stepsize,
    stepsize.decayrate = 0.6,
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
  params$Theta.avg <- params$Theta.avg[,attributes(PhiAvgEval)$match_id,,,drop = F]
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
      Method = rep(paste0("Pspline-", optMethod), nRecord + 1),
      StepSize = rep(stepsize, nRecord + 1),
      N = c(0, nIters),
      nBatch = rep(nBatch, nRecord + 1),
      Time = times,
      # AV = c(NA,vscores[,"av"]),
      # EWMABV = c(NA,vscores[,"ewmabv"]),
      RMSEphi1 = rmseAll[, 1],
      RMSEphi2 = rmseAll[, 2],
      RMSEphi3 = rmseAll[, 3],
      RMSEphi1.avg = rmseAll.avg[, 1],
      RMSEphi2.avg = rmseAll.avg[, 2],
      RMSEphi3.avg = rmseAll.avg[, 3]
    )
  )

  for (ip in seq_len(nPass)) {
    ii <- ip * nRecord.1pass
    ThetaRecord[[ip]][[optMethod]] <- ThetaAll.avg[,,ii + 1]
  }
}


## Batch FPCA ------------------------
message("Start Batch FPCA: RegFPCA ---------------")
batch_start <- Sys.time()
fitBatch <- fpca.reg(
  dat$Ly,
  dat$Ltid,
  inits = inits,
  meanfun = FALSE,
  npc = npc,
  maxIter = 100,
  nu = 0.5,
  verbose = FALSE,
  record_iterations = TRUE,
  tau = 10^seq(-9, -4, 1),
  use_validation_set = FALSE,
  refine_alpha = FALSE
)
ThetaBatch <- fitBatch$Theta
muBatchEval <- eval_fd(evalGrid, FuncData(fitBatch$theta_mu, basis))
PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
PhiBatchEval <- match_fpc(PhiBatchEval, PhiTrueEval[, 1:q, drop = F])
batch_end <- Sys.time()
rmseBatch <- sqrt(colMeans((PhiBatchEval - PhiTrueEval[, 1:q])^2))
res <- rbind(
  res,
  data.frame(
    seed = seed,
    Method = "Batch-SOAP",
    StepSize = NA,
    N = N,
    nBatch = N,
    Time = difftime(batch_end, batch_start, units = 'secs') +
      difftime(end_time.init, start_time.init, units = 'secs'),
    # AV = fit0$av,
    # EWMABV = fit0$ewmabv,
    RMSEphi1 = rmseBatch[1],
    RMSEphi2 = rmseBatch[2],
    RMSEphi3 = rmseBatch[3],
    RMSEphi1.avg = rmseBatch[1],
    RMSEphi2.avg = rmseBatch[2],
    RMSEphi3.avg = rmseBatch[3]
  )
)

# End experiment ----------------------------------------------------------

readr::write_csv(res, file.path(dirpath, filename))
saveRDS(ThetaRecord, file=file.path(dirpath, paste0("Theta_",exprmt,"_sd",seedtext,".rds")))

message("Finishing replication: ", seed)
set.seed(NULL)
