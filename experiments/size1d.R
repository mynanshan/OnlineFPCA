#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=size1d
#SBATCH --output=logs/size1d_%j.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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
parser$add_argument(
  "--ninit",
  type = "integer",
  default = 200,
  help = "Number of subjects for initialization."
)
parser$add_argument(
  "--algo",
  type = "integer",
  default = 1,
  choices = c(1, 2, 3, 4, 5),
  help = "Test what method (1: OnlineFPCA, 2: PACE, 3: FACE, 4: REML, 5: SOAP)"
)
args <- parser$parse_args()
seed <- as.integer(args[["seed"]])
N <- as.integer(args[["N"]])
Ninit <- as.integer(args[["ninit"]])
algo <- as.integer(args[["algo"]])

model_name <- c(
  "1" = "OnlineFPCA",
  "2" = "PACE",
  "3" = "FACE",
  "4" = "REML",
  "5" = "SOAP"
)[as.character(algo)]

cat("========= Start Experiment ( Seed =", seed, ")\n")
cat("Algo: ", algo, "\n")
cat("Model: ", model_name, "\n")

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

noiseList <- c(0.5)
nBatch <- 5
nParams <- 6
nPass <- 1
nBlock <- 100
initMethod <- "face"
N <- 25000
Nsmalls <- c(200, 500, seq(1000, N, 1000))
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsizeList <- c(0.1)
stepsize.min <- 1e-3
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune

exprmt <- "size1d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed <- ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("simu_", exprmt, "_alg", algo, "_sd", seedtext, ".csv")

res <- data.frame(
  noise = numeric(),
  seed = numeric(),
  Method = character(),
  StepSize = numeric(),
  N = numeric(),
  Ninit = numeric(),
  initMethod = character(),
  npc = numeric(),
  nBatch = numeric(),
  Time = numeric(),
  RMSEphi1 = numeric(),
  RMSEphi2 = numeric(),
  RMSEphi3 = numeric(),
  RMSEphi1.avg = numeric(),
  RMSEphi2.avg = numeric(),
  RMSEphi3.avg = numeric(),
  tau = numeric()
)

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
  
  message("Start initialization.")

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
  g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init) / q)

  inits <- list(
    Theta = ThetaInit,
    lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
    sigma2 = sigma2Init
  )
  
  if (algo == 1) {
    # Online FPCA  ------------------------------------------------------

    message("Starting model: OnlineFPCA.")
  
    for (stepsize0 in stepsizeList) {
      sgd_lr0 <- stepsize0 / g_Theta_init_norm
      message("> Step size = ", stepsize0)
  
      for (sgdtype in c("sgd", "adagrad")) {
        message(">>> Using model: OnlineFPCA-", sgdtype)
  
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
          adareset = 25 * nBlockIter,
          adareset.end = 50 * nBlockIter,
          asgd.start = 10 * nBlockIter + 1,
          asgd.reset = 20 * nBlockIter,
          asgd.reset.end = 50 * nBlockIter,
          asgd.end = Inf,
          nIter.1stTune = nRoundNoTune * nBlockIter,
          nIter.lastTune = nIter1pass,
          nIter.tauNoIncrease = floor(0.3 * nIter1pass),
          period.tune = nBlockIter,
          period.record = nBlockIter,
          verbose = FALSE
        )
  
        check <- fit$check

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
          drop = FALSE
        ]
        params$Theta[, attributes(PhiEstEval)$flipped, , ] <-
          -params$Theta[, attributes(PhiEstEval)$flipped, , ]
        PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg, basis))
        PhiAvgEval <- match_fpc(PhiAvgEval, PhiTrueEval[, 1:q, drop = F])
        params$Theta.avg <- params$Theta.avg[,
          attributes(PhiAvgEval)$match_id, , ,
          drop = FALSE
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

        times <- if (nPass == 1) {
          colSums(fit$time.history[, 1:nIter1pass])
        } else {
          c(
            colSums(fit$time.history[, 1:nIter1pass]),
            fit$time.history[1, (nIter1pass + 1):(nPass * nIter1pass)]
          )
        }
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
            Ninit = rep(Ninit, nRecord + 1),
            initMethod = rep(initMethod, nRecord + 1),
            npc = rep(q, nRecord + 1),
            nBatch = rep(nBatch, nRecord + 1),
            Time = times,
            RMSEphi1 = rmseAll[, 1],
            RMSEphi2 = rmseAll[, 2],
            RMSEphi3 = rmseAll[, 3],
            RMSEphi1.avg = rmseAll.avg[, 1],
            RMSEphi2.avg = rmseAll.avg[, 2],
            RMSEphi3.avg = rmseAll.avg[, 3],
            tau = tau_path$tau_path[seq_len(rmseAll)]
          )
        )
  
      } # end of all sgdtype
    } # end of all step sizes
  } # endif algo==1

  if (algo == 2) {
    message("Starting model: PACE.")
    for (Nsmall in Nsmalls) {
      # Batch method 1: PACE ----------------------------------------------------
      message("Using model: PACE; Nsmall = ", Nsmall)
      batch_start <- Sys.time()
      fitBatch <- fdapace::FPCA(
        dat$Ly[1:Nsmall],
        dat$Lt[1:Nsmall],
        optns = list(
          dataType = "Sparse",
          nRegGrid = neval,
          userMu = list(t = evalGrid, mu = muTrueEval),
          methodBwCov = "GMeanAndGCV"
        )
      )
      muBatchEval <- fitBatch$mu
      theta_muBatch <- smooth_basis(
        evalGrid,
        muBatchEval,
        basis,
        lambda = 1e-10
      )$coefs
      eig <- RSpectra::eigs_sym(fitBatch$smoothedCov, npc)
      PhiBatchEval <- flip_direc(eig$vectors, PhiTrueEval[, 1:npc])
      ThetaBatch <- smooth_basis(
        evalGrid,
        PhiBatchEval,
        basis,
        lambda = 1e-10
      )$coefs
      norm_factor <- diag(t(ThetaBatch) %*% G %*% ThetaBatch)
      ThetaBatch <- sweep(ThetaBatch, 2, sqrt(norm_factor), "/")
      lambdaBatch <- eig$values
      sigma2Batch <- fitBatch$sigma2
      batch_end <- Sys.time()
      PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
      PhiBatchEval <- match_fpc(PhiBatchEval, PhiTrueEval[, 1:npc, drop = F])
      rmseBatch <- sqrt(colMeans((PhiBatchEval - PhiTrueEval[, 1:npc])^2))
      res <- rbind(
        res,
        data.frame(
          noise = noise_sd,
          seed = seed,
          Method = "Batch-PACE",
          StepSize = NA,
          N = Nsmall,
          Ninit = NA,
          initMethod = NA,
          npc = npc,
          nBatch = Nsmall,
          Time = difftime(batch_end, batch_start, units = "secs"),
          RMSEphi1 = rmseBatch[1],
          RMSEphi2 = rmseBatch[2],
          RMSEphi3 = rmseBatch[3],
          RMSEphi1.avg = rmseBatch[1],
          RMSEphi2.avg = rmseBatch[2],
          RMSEphi3.avg = rmseBatch[3],
          tau = NA
        )
      )
    }
  } # endif algo==2

  if (algo == 3) {
    message("Starting model: FACE.")
    for (Nsmall in Nsmalls) {
      # Batch method 2: FACE ----------------------------------------------------
      message("Using model: FACE; Nsmall = ", Nsmall)
      tmpdat <- data.frame(
        argvals = unlist(dat$Lt[1:Nsmall]),
        subj = rep(1:Nsmall, dat$Lmi[1:Nsmall]),
        y = unlist(dat$Ly[1:Nsmall])
      )
      batch_start <- Sys.time()
      fitBatch <- face::face.sparse(
        tmpdat,
        center = FALSE,
        argvals.new = evalGrid,
        knots = p
      )
      muBatchEval <- fitBatch$mu.new
      theta_muBatch <- smooth_basis(
        evalGrid,
        muBatchEval,
        basis,
        lambda = 1e-10
      )$coefs
      eig <- RSpectra::eigs_sym(fitBatch$Chat.new, npc)
      PhiBatchEval <- flip_direc(eig$vectors, PhiTrueEval[, 1:npc])
      ThetaBatch <- smooth_basis(
        evalGrid,
        PhiBatchEval,
        basis,
        lambda = 1e-10
      )$coefs
      lambdaBatch <- eig$values
      sigma2Batch <- fitBatch$sigma2
      norm_factor <- diag(t(ThetaBatch) %*% G %*% ThetaBatch)
      ThetaBatch <- sweep(ThetaBatch, 2, sqrt(norm_factor), "/")
      lambdaBatch <- lambdaBatch * norm_factor
      PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
      batch_end <- Sys.time()
      PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
      PhiBatchEval <- match_fpc(PhiBatchEval, PhiTrueEval[, 1:npc, drop = F])
      rmseBatch <- sqrt(colMeans((PhiBatchEval - PhiTrueEval[, 1:npc])^2))
      res <- rbind(
        res,
        data.frame(
          noise = noise_sd,
          seed = seed,
          Method = "Batch-FACE",
          StepSize = NA,
          N = Nsmall,
          Ninit = NA,
          initMethod = NA,
          npc = npc,
          nBatch = Nsmall,
          Time = difftime(batch_end, batch_start, units = "secs"),
          # AV = fit0$av,
          # EWMABV = fit0$ewmabv,
          RMSEphi1 = rmseBatch[1],
          RMSEphi2 = rmseBatch[2],
          RMSEphi3 = rmseBatch[3],
          RMSEphi1.avg = rmseBatch[1],
          RMSEphi2.avg = rmseBatch[2],
          RMSEphi3.avg = rmseBatch[3],
          tau = NA
        )
      )
    }
  } # endif algo==3

  if (algo == 4) {
    message("Starting model: REML.")
    for (Nsmall in Nsmalls) {
      # Batch method 3: REML ----------------------------------------------------
      message("Using model: REML; Nsmall = ", Nsmall)
      tmpdat <- cbind(
        rep(1:Nsmall, dat$Lmi[1:Nsmall]),
        unlist(dat$Ly[1:Nsmall]),
        unlist(dat$Lt[1:Nsmall])
      )
      batch_start <- Sys.time()
      fitBatch <- fpca::fpca.mle(
        tmpdat,
        M.set = c(5:9),
        r.set = npc,
        ini.method = "loc",
        basis.method = "bs",
        max.step = 50,
        grid.l = evalGrid,
        grids = evalGrid
      )
      muBatchEval <- fitBatch$fitted_mean
      theta_muBatch <- smooth_basis(
        evalGrid,
        muBatchEval,
        basis,
        lambda = 1e-10
      )$coefs
      PhiBatchEval <- flip_direc(t(fitBatch$eigenfunctions), PhiTrueEval[, 1:npc])
      ThetaBatch <- smooth_basis(
        evalGrid,
        PhiBatchEval,
        basis,
        lambda = 1e-10
      )$coefs
      lambdaBatch <- fitBatch$eigenvalues
      norm_factor <- diag(t(ThetaBatch) %*% G %*% ThetaBatch)
      ThetaBatch <- sweep(ThetaBatch, 2, sqrt(norm_factor), "/")
      lambdaBatch <- lambdaBatch * norm_factor
      PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
      sigma2Batch <- unname(fitBatch$error_var)
      batch_end <- Sys.time()
      PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
      PhiBatchEval <- match_fpc(PhiBatchEval, PhiTrueEval[, 1:npc, drop = F])
      rmseBatch <- sqrt(colMeans((PhiBatchEval - PhiTrueEval[, 1:npc])^2))
      res <- rbind(
        res,
        data.frame(
          noise = noise_sd,
          seed = seed,
          Method = "Batch-REML",
          StepSize = NA,
          N = Nsmall,
          Ninit = NA,
          initMethod = NA,
          npc = npc,
          nBatch = Nsmall,
          Time = difftime(batch_end, batch_start, units = "secs"),
          # AV = fit0$av,
          # EWMABV = fit0$ewmabv,
          RMSEphi1 = rmseBatch[1],
          RMSEphi2 = rmseBatch[2],
          RMSEphi3 = rmseBatch[3],
          RMSEphi1.avg = rmseBatch[1],
          RMSEphi2.avg = rmseBatch[2],
          RMSEphi3.avg = rmseBatch[3],
          tau = NA
        )
      )
    }
  } # endif algo==4

  if (algo == 5) {
    message("Starting model: SOAP.")
    inits_soap <- inits
    for (Nsmall in Nsmalls) {
      # Batch method 4: SOAP ------------------------------------------
      message("Using model: SOAP; Nsmall = ", Nsmall)
      batch_start <- Sys.time()
      fitBatch <- fpca.reg(
        dat$Ly[1:Nsmall],
        dat$Ltid[1:Nsmall],
        inits_soap,
        meanfun = FALSE,
        npc = 3,
        maxIter = 100,
        nu = 1,
        verbose = FALSE,
        record_iterations = TRUE,
        tau = 10^seq(-3, -1, 0.5),
        use_validation_set = FALSE,
        refine_alpha = FALSE
      )
      ThetaBatch <- fitBatch$Theta
      muBatchEval <- eval_fd(evalGrid, FuncData(fitBatch$theta_mu, basis))
      PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
      PhiBatchEval <- match_fpc(PhiBatchEval, PhiTrueEval[, 1:npc, drop = F])
      batch_end <- Sys.time()
      rmseBatch <- sqrt(colMeans((PhiBatchEval - PhiTrueEval[, 1:npc])^2))
      res <- rbind(
        res,
        data.frame(
          noise = noise_sd,
          seed = seed,
          Method = "Batch-SOAP",
          StepSize = NA,
          N = Nsmall,
          Ninit = NA,
          initMethod = NA,
          npc = npc,
          nBatch = Nsmall,
          Time = difftime(batch_end, batch_start, units = "secs") +
            difftime(end_time.init, start_time.init, units = "secs"),
          # AV = fit0$av,
          # EWMABV = fit0$ewmabv,
          RMSEphi1 = rmseBatch[1],
          RMSEphi2 = rmseBatch[2],
          RMSEphi3 = rmseBatch[3],
          RMSEphi1.avg = rmseBatch[1],
          RMSEphi2.avg = rmseBatch[2],
          RMSEphi3.avg = rmseBatch[3],
          tau = NA
        )
      )
      inits_soap$Theta <- ThetaBatch
      inits_soap$alphamat <- fitBatch$alphamat
    }
  } # endif algo==5
}


# End experiment ----------------------------------------------------------

readr::write_csv(res, file.path(dirpath, filename))

message("Finishing replication: ", seed)
set.seed(NULL)
