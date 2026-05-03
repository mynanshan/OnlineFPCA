#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=size2d
#SBATCH --output=logs/size2d_%j.out
#SBATCH --time=12:00:00
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
parser$add_argument(
  "--compare",
  type = "integer",
  default = 1,
  choices = c(0, 1),
  help = "Whether to run competing methods."
)
parser$add_argument(
  "--simple",
  type = "integer",
  default = 0,
  choices = c(0, 1),
  help = "Whether to simplify experiment settings."
)
args <- parser$parse_args()
seed <- as.integer(args[["seed"]])
compare <- as.logical(args[["compare"]])
simple <- as.logical(args[["simple"]])

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

noiseList <- c(0.1, 0.5, 1.)
nBatch <- 5
nParams <- 6
nPass <- 5
nBlock <- 100
Ninit <- 100
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
# stepsizeList <- c(0.2, 0.1, 0.05)
stepsizeList <- c(0.1)
# stepsize0 <- 0.1
stepsize.min <- 1e-4
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune

exprmt <- "size2d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed <- ceiling(log10(seed + 1))
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

G <- get_basis_inprod_matrix(basis)
GR <- chol(G)
invG <- chol2inv(GR)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv

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

  # Online FPCA  ------------------------------------------------------------

  message("Start OnlineFPCA ---------------")

  start_time.init <- Sys.time()
  Yinit <- dat$Ly[1:Ninit]
  Tinit <- dat$Lt[1:Ninit]
  LmiInit <- dat$Lmi[1:Ninit]
  TidInit <- dat$Ltid[1:Ninit]
  nsub <- round((LmiInit + runif(Ninit, min = -0.1, max = 0.1)) * 0.5)
  nsub <- pmax(nsub, 4)
  nsub <- pmin(nsub, 20)
  nsub <- pmin(nsub, LmiInit)
  sampleIds <- mapply(
    \(mi, mi_sub) {
      sort(sample(1:mi, mi_sub))
    },
    LmiInit,
    nsub
  )
  Yinit <- unlist(mapply(
    \(yi, ids) {
      yi[ids]
    },
    Yinit,
    sampleIds
  ))
  Tinit <- do.call(
    rbind,
    mapply(
      \(ti, ids) {
        ti[ids, ]
      },
      Tinit,
      sampleIds
    )
  )
  TidInit <- unlist(mapply(
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
  FpcaOut <- fpca.mOpCov(OUT = CovRes)
  lambdaInit <- FpcaOut$Eigen$values[1:q]
  PhiInitEvalFull <- computeEigen(evalGrid, CovRes, FpcaOut)
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
  innerIds <- which(
    evalGrid[, 1] > 0.1 &
      evalGrid[, 1] < 0.9 &
      evalGrid[, 2] > 0.1 &
      evalGrid[, 2] < 0.9
  )
  innerPoints <- which(TidInit %in% innerIds)
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
  g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init) / q)

  inits <- list(
    Theta = ThetaInit,
    lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
    sigma2 = sigma2Init
  )

  for (stepsize0 in stepsizeList) {
    sgd_lr0 <- stepsize0 / g_Theta_init_norm
    message("> Step size = ", stepsize0)

    for (sgdtype in c("sgd", "adagrad")) {
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
        verbose = TRUE
      )

      tau.min <- fit$tau.min
      l <- which(fit$tau == tau.min)[1]
      Theta <- fit$Theta[, , l]
      lambda <- fit$lambda[, l]
      sigma2 <- fit$sigma2[l]
      Theta.avg <- fit$Theta.avg[, , l]
      lambda.avg <- fit$lambda.avg[, l]
      sigma2.avg <- fit$sigma2.avg[l]

      params <- fit$params.history$params
      vcrits <- fit$vcrit.history

      # check eigenfunctions
      PhiEstEval <- eval_fd(evalGrid, FuncData(Theta, basis))
      PhiEstEval <- match_fpc(PhiEstEval, PhiTrueEval[, 1:q, drop = F])
      params$Theta <- params$Theta[,
        attributes(PhiEstEval)$match_id, , ,
        drop = F
      ]
      params$Theta[, attributes(PhiEstEval)$flipped, , ] <-
        -params$Theta[, attributes(PhiEstEval)$flipped, , ]
      PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg, basis))
      PhiAvgEval <- match_fpc(PhiAvgEval, PhiTrueEval[, 1:q, drop = F])
      params$Theta.avg <- params$Theta.avg[,
        attributes(PhiAvgEval)$match_id, , ,
        drop = F
      ]
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
        \(i) params$Theta[, , tau_path_id_extend[i], i],
        simplify = "array"
      )
      rmseAll <- rmse_phi(ThetaAll, phiTrueFunc$coefs, B)
      ThetaAll.avg <- sapply(
        seq_along(fit$params.history$iter.params),
        \(i) params$Theta.avg[, , tau_path_id_extend[i], i],
        simplify = "array"
      )
      rmseAll.avg <- rmse_phi(ThetaAll.avg, phiTrueFunc$coefs, B)

      times <- c(
        colSums(fit$time.history[, 1:nIter1pass]),
        fit$time.history[1, (nIter1pass + 1):(nPass * nIter1pass)]
      )
      times <- colSums(matrix(times, nrow = nBlockIter))
      times <- c(difftime(end_time.init, start_time.init, units = "secs"), times)

      res <- rbind(
        res,
        data.frame(
          noise = noise_sd,
          seed = rep(seed, nRecord + 1),
          Method = rep(paste0("OnlineFPCA-", sgdtype), nRecord + 1),
          StepSize = rep(stepsize0, nRecord + 1),
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

      if (seed == 1234 && noise_sd == 0.1) {
        pdf(
          file.path(dirpath, paste0("taupath-sim2d-", sgdtype, ".pdf")),
          width = 8, height = 5
        )
        do.call(
          plot.tau_path2,
          list(
            tau.history = fit$tau.select$tau.history,
            tau.selectId = fit$tau.select$tau.selectId,
            final.id = l, nbatch_per_block = nBlock / nBatch
          )
        )
        dev.off()
      }
    } # end of all sgdtype
  } # end of all step sizes

  ## Batch FPCA ------------------------
  if (compare) {

    # SOAP -------------------
    batch_start <- Sys.time()
    fitBatch <- fpca.reg(
      dat$Ly,
      dat$Ltid,
      inits = list(
        Theta = ThetaInit,
        lambda = lambdaInit,
        sigma2 = sigma2Init
      ),
      meanfun = FALSE,
      npc = q,
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
        noise = noise_sd,
        seed = seed,
        Method = "Batch-SOAP",
        StepSize = NA,
        N = N,
        nBatch = N,
        Time = difftime(batch_end, batch_start, units = "secs") +
          difftime(end_time.init, start_time.init, units = "secs"),
        RMSEphi1 = rmseBatch[1],
        RMSEphi2 = rmseBatch[2],
        RMSEphi3 = rmseBatch[3],
        RMSEphi1.avg = rmseBatch[1],
        RMSEphi2.avg = rmseBatch[2],
        RMSEphi3.avg = rmseBatch[3]
      )
    )

    # mOpCov -----------------------
    # Nsmalls <- c(200, 300, 400, 500)
    Nsmalls <- c(200, 500)
    for (Nsmall in Nsmalls) {
      batch_start <- Sys.time()
      CovRes <- mOpCov(
        location = do.call(rbind, dat$Lt[1:Nsmall]),
        x = unlist(dat$Ly[1:Nsmall]),
        subject = rep(1:Nsmall, times = dat$Lmi[1:Nsmall]),
        q = c(6, 6),
        lam = list(lam = 1e-10, alpha = 1e-6),
        ker = "cos"
      )
      FpcaOut <- fpca.mOpCov(OUT = CovRes)
      lambdaBatch <- FpcaOut$Eigen$values[1:q]
      PhiBatchEvalFull <- computeEigen(evalGrid, CovRes, FpcaOut)
      PhiBatchEval <- PhiBatchEvalFull[, 1:q]
      ThetaBatch <- smooth_basis(
        evalGrid,
        PhiBatchEval,
        basis,
        lambda = 1e-8
      )$coefs
      norm_factor <- diag(t(ThetaBatch) %*% G %*% ThetaBatch)
      ThetaBatch <- sweep(ThetaBatch, 2, sqrt(norm_factor), "/")
      lambdaBatch <- lambdaBatch * norm_factor
      PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
      PhiBatchEval <- flip_direc(PhiBatchEval, PhiTrueEval[, 1:q])
      flipId <- attributes(PhiBatchEval)$flipped
      ThetaBatch[, flipId] <- -ThetaBatch[, flipId]
      batch_end <- Sys.time()
      
      rmseBatch <- sqrt(colMeans((PhiBatchEval - PhiTrueEval[, 1:q])^2))
      res <- rbind(
        res,
        data.frame(
          noise = noise_sd,
          seed = seed,
          Method = paste0("Batch-mOpCov-N", Nsmall),
          StepSize = NA,
          N = Nsmall,
          nBatch = N,
          Time = difftime(batch_end, batch_start, units = "secs"),
          RMSEphi1 = rmseBatch[1],
          RMSEphi2 = rmseBatch[2],
          RMSEphi3 = rmseBatch[3],
          RMSEphi1.avg = rmseBatch[1],
          RMSEphi2.avg = rmseBatch[2],
          RMSEphi3.avg = rmseBatch[3]
        )
      )
    } # end all Nsmall
  } # end competing method
}


# End experiment ----------------------------------------------------------

if (!simple) readr::write_csv(res, file.path(dirpath, filename))

message("Finishing replication: ", seed)
set.seed(NULL)
