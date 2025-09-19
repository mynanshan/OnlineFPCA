#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=boot-fpca2d
#SBATCH --output=logs/boot-fpca2d_%j.out
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
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

mOpCov_path <- "external_codes/mOpCov/"
source(paste0(mOpCov_path, "mOpCov_prep.R"))
Rcpp::sourceCpp(paste0(mOpCov_path, "mOpCov_cpp.cpp"))

nBatch <- 5
noise_sd <- 0.5
alpha <- 2
nParams <- 6
nPass <- 3
nBlock <- 100
Ninit <- 100
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize <- 2e-1
stepsize.min <- 1e-1
sgd.step.scale <- 1 # use a smaller step size for sgd
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
asgd.use <- TRUE

exprmt <- "boot2d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
basename <- paste0(exprmt, "_sd", seedtext)

res <- data.frame(
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

tau_paths_mat <- c()

rp <- get_rp.wang2020(alpha = alpha)
m_min <- NULL
m_max <- NULL
m_mean <- 25
m_sd <- 6
m_type <- "gaussian"

t0 <- rp$t0
t1 <- rp$t1
D <- rp$ndim

q <- 3

# evaluate true mean and eigenfunctions
nevalList <- c(51,51)
neval <- prod(nevalList)
evalGridList <- lapply(1:rp$ndim, \(d) seq(t0[d], t1[d], length.out = nevalList[d]))
evalGrid <- margins2grid(evalGridList)
muTrueEval <- rp$meanfun(evalGrid)
PhiTrueEval <- rp$eigfun(evalGrid)
lambdaTrue <- rp$eigval
covTrueEval <- PhiTrueEval %*% diag(lambdaTrue) %*% t(PhiTrueEval)
evalGridListSmall <- lapply(1:rp$ndim, \(d) seq(t0[d], t1[d], length.out = 21))
evalGridSmall <- margins2grid(evalGridListSmall)

set.seed(1234)

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

set.seed(seed)
ord = sample(1:N, N, replace = TRUE)
dat$Ly <- dat$Ly[ord]
dat$Lt <- dat$Lt[ord]
dat$Ltid <- dat$Ltid[ord]
dat$Lmi <- dat$Lmi[ord]
dat$score <- dat$score[ord,]

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
nbasis1d <- 7
basis <- TensorBasis(list(
  create.bspline.basis(c(t0[1],t1[1]), nbasis = nbasis1d, norder = 4),
  create.bspline.basis(c(t0[2],t1[2]), nbasis = nbasis1d, norder = 4)))
p <- attr(basis, "nbasis")

iRecord <- c(10, 30, 50, 100, 150)
ThetaRecord <- rep(list(list(SGD = NULL, Adam = NULL)), length(iRecord))

muTrueFunc <- smooth_basis(evalGrid, muTrueEval, basis)
phiTrueFunc <- smooth_basis(evalGrid, PhiTrueEval, basis)
rmseMuBest <- Metrics::rmse(muTrueEval, eval_fd(evalGrid, muTrueFunc))
rmsePhiBest <- Metrics::rmse(PhiTrueEval, eval_fd(evalGrid, phiTrueFunc))

B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)


# Online FPCA  ------------------------------------------------------------

message("Start OnlineFPCA ---------------")

start_time.init <- Sys.time()
Yinit = dat$Ly[1:Ninit]
Tinit = dat$Lt[1:Ninit]
LmiInit = dat$Lmi[1:Ninit]
TidInit <- dat$Ltid[1:Ninit]
# do a subsampling within each subject to accelerate mOpCov
nsub = round((LmiInit + runif(Ninit,min=-0.1,max=0.1)) * 0.5)
nsub = pmax(nsub, 4)
nsub = pmin(nsub, 20)
nsub = pmin(nsub, LmiInit)
sampleIds = mapply(\(mi, mi_sub) {sort(sample(1:mi, mi_sub))}, LmiInit, nsub)
Yinit = unlist(mapply(\(yi, ids) {yi[ids]}, Yinit, sampleIds))
Tinit = do.call(rbind, mapply(\(ti, ids) {ti[ids,]}, Tinit, sampleIds))
TidInit = unlist(mapply(\(tid, ids) {tid[ids]}, TidInit, sampleIds))
subject_ids <- rep(1:Ninit, times = nsub)
CovRes <- mOpCov(location=Tinit, x=Yinit, subject=subject_ids,
                  q=c(6,6), lam = list(lam = 1e-10, alpha = 1e-6), ker = "cos")
FpcaOut <- fpca.mOpCov(OUT = CovRes)
lambdaInit <- FpcaOut$Eigen$values[1:q]
PhiInitEvalFull <- computeEigen(evalGrid,CovRes,FpcaOut)
PhiInitEval <- PhiInitEvalFull[,1:q]
ThetaInit <- smooth_basis(evalGrid, PhiInitEval, basis, lambda = 1e-8)$coefs
norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- eval_fd(evalGrid, FuncData(ThetaInit, basis))
PhiInitEval <- flip_direc(PhiInitEval, PhiTrueEval[,1:q])
flipId <- attributes(PhiInitEval)$flipped
ThetaInit[,flipId] <- -ThetaInit[,flipId]
# fit diagonal element, estimate sigma2
innerIds = which(evalGrid[,1] > 0.1 & evalGrid[,1] < 0.9 &
                    evalGrid[,2] > 0.1 & evalGrid[,2] < 0.9)
innerPoints = which(TidInit %in% innerIds)
covInitEvalDiag <- as.vector(PhiInitEvalFull^2 %*% FpcaOut$Eigen$values)
sigma2Init <- max(mean((Yinit^2 - covInitEvalDiag[TidInit])[innerPoints]), 1e-3)
end_time.init <- Sys.time()

nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)

inits <- list(Theta = ThetaInit, lambda = lambdaInit, sigma2 = sigma2Init)

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
      maxtau = 1e-3,
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
      Method = rep(paste0("Pspline-", optMethod), nRecord + 1),
      StepSize = rep(stepsize, nRecord + 1),
      N = c(0, nIters),
      nBatch = rep(nBatch, nRecord + 1),
      Time = times,
      RMSEphi1 = rmseAll[, 1],
      RMSEphi2 = rmseAll[, 2],
      RMSEphi3 = rmseAll[, 3],
      RMSEphi1.avg = rmseAll.avg[, 1],
      RMSEphi2.avg = rmseAll.avg[, 2],
      RMSEphi3.avg = rmseAll.avg[, 3]
    )
  )

  tau_paths_mat <- cbind(tau_paths_mat, tau_path)
  colnames(tau_paths_mat)[ncol(tau_paths_mat)] <- optMethod

  for (ir in seq_along(iRecord)) {
    ii <- iRecord[ir]
    ThetaRecord[[ir]][[optMethod]] <- ThetaAll.avg[,,ii + 1]
  }
}

# End experiment ----------------------------------------------------------

readr::write_csv(res, file.path(dirpath, paste0("simu_", basename, ".csv")))
saveRDS(tau_paths_mat, file.path(dirpath, paste0("taupath_", basename, ".csv")))
saveRDS(ThetaRecord, file=file.path(dirpath, paste0("Theta_",exprmt,"_sd",seedtext,".rds")))

message("Finishing replication: ", seed)
set.seed(NULL)
