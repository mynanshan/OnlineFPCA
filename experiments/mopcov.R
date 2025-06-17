arg <- commandArgs(trailingOnly = TRUE)
seed0 <- as.numeric(arg[1])
if (is.na(seed0)) stop("UNKNOWN SEED")

set.seed(seed0)
# set.seed(1)


mOpCov_path <- "external_codes/mOpCov/"
source(paste0(mOpCov_path, "mOpCov_prep.R"))
Rcpp::sourceCpp(paste0(mOpCov_path, "mOpCov_cpp.cpp"))

library(ManifoldOptim)

## Experiment settings
nRep <- 20
# snrList <- c(2,5)
# alphaList <- c(1.2,2,3)
# snrList <- c(2)
alphaList <- c(2)
# noiseList <- c(0.1, 0.4)
noiseList <- c(0.1)
N <- 5000
NbatchList = c(1000, 2000, 5000)

dirpath <- "experiments"
exprmt <- paste0("mopcov", "N", N)

if (!dir.exists(dirpath)) dir.create(dirpath)

#%%%%%%%%%%%%%%%%%%%%%%%%%


library(fda)
library(Matrix)


source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./data_generation/generator.R")

# set.seed(2501)

for (alpha in alphaList) {
  for (noise_sd in noiseList) {
    message("alpha=", alpha, ", ", "noise_sd=", noise_sd)
    
    filename <- file.path(
      dirpath, paste0(exprmt, "_alpha", stringr::str_replace(alpha, "\\.", "p"),
                      "_noise", stringr::str_replace(noise_sd, "\\.", "p"), "_",
                      stringr::str_remove_all(
                        as.character(round(as.numeric(Sys.time()))), "[ \\-:]"), ".csv"))
    
    res <- data.frame(iRep = numeric(),
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
                      RMSEphi3.avg = numeric())
    
    rp <- get_rp.wang2020(alpha = alpha)
    # m_min <- 5; m_max <- 15; m_mean <- NULL; m_sd <- NULL
    # m_type <- "unif"
    m_min <- NULL; m_max <- NULL; m_mean <- 25; m_sd <- 6
    m_type <- "gaussian"
    
    t0 <- rp$t0; t1 <- rp$t1
    D <- rp$ndim
    
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
    
    q <- 3 # number of PCs to be extracted
    
    for (ir in 1:nRep) {
      message("iter ", ir)
      
      dat <- get_measurements(
        rp, n = N, m_min = m_min, m_max = m_max, m_mean = m_mean, m_sd = m_sd,
        design_type = "random", m_type = m_type, sigma = noise_sd)
      fdata_generator <- function(n, total_count) {
        idx <- (total_count):(total_count+n-1) %% N + 1
        return(list(Ly=dat$Ly[idx], Ltid=dat$Ltid[idx], Lt=dat$Lt[idx], Lmi=dat$Lmi[idx]))
      }
      tgrid <- dat$tgrid
      
      # set B-spline basis
      nbasis1d <- 7
      basis <- TensorBasis(list(
        create.bspline.basis(c(t0[1],t1[1]), nbasis = nbasis1d, norder = 4),
        create.bspline.basis(c(t0[2],t1[2]), nbasis = nbasis1d, norder = 4)))
      p <- attr(basis, "nbasis")
      
      # best possible basis approximations
      muTrueFunc <- smooth_basis(evalGrid, muTrueEval, basis)
      phiTrueFunc <- smooth_basis(evalGrid, PhiTrueEval, basis)
      rmseMuBest <- Metrics::rmse(muTrueEval, eval_fd(evalGrid, muTrueFunc))
      rmsePhiBest <- Metrics::rmse(PhiTrueEval, eval_fd(evalGrid, phiTrueFunc))
      # Metrics::rmse(PhiTrueEval[,1:q], eval_fd(evalGrid, phiTrueFunc)[,1:q])
      
      for (Nbatch in NbatchList) {
        # mOpCov --------------------------------
        batch_start <- Sys.time()
        Ybatch = unlist(dat$Ly[1:Nbatch])
        Tbatch = do.call(rbind, dat$Lt[1:Nbatch])
        LmiBatch = unlist(dat$Lmi[1:Nbatch])
        # TidBatch <- (dat$Ltid[1:Nbatch])
        subject_ids <- rep(1:Nbatch, times = LmiBatch)
        CovRes <- mOpCov(location=Tbatch, x=Ybatch, subject=subject_ids,
                         q=c(6,6), lam = list(lam = 1e-10, alpha = 1e-6), ker = "cos")
        FpcaOut <- fpca.mOpCov(OUT = CovRes)
        batch_end <- Sys.time()
        lambdaBatch <- FpcaOut$Eigen$values[1:q]
        PhiBatchEvalFull <- computeEigen(evalGrid,CovRes,FpcaOut)
        PhiBatchEval <- PhiBatchEvalFull[,1:q]
        ThetaBatch <- smooth_basis(evalGrid, PhiBatchEval, basis, lambda = 1e-8)$coefs
        norm_factor <- diag(t(ThetaBatch) %*% G %*% ThetaBatch)
        ThetaBatch <- sweep(ThetaBatch, 2, sqrt(norm_factor), "/")
        lambdaBatch <- lambdaBatch * norm_factor
        PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
        PhiBatchEval <- flip_direc(PhiBatchEval, PhiTrueEval[,1:q])
        flipId <- attributes(PhiBatchEval)$flipped
        ThetaBatch[,flipId] <- -ThetaBatch[,flipId]
        # fit diagonal element, estimate sigma2
        # innerIds = which(evalGrid[,1] > 0.1 & evalGrid[,1] < 0.9 &
        #                    evalGrid[,2] > 0.1 & evalGrid[,2] < 0.9)
        # innerPoints = which(TidBatch %in% innerIds)
        # covBatchEvalDiag <- as.vector(PhiBatchEvalFull^2 %*% FpcaOut$Eigen$values)
        # sigma2Batch <- mean((Ybatch^2 - covBatchEvalDiag[TidBatch])[innerPoints])
        rmseBatch <- sqrt(colMeans((PhiBatchEval - PhiTrueEval[,1:q])^2))
        res <- rbind(res, 
                     data.frame(iRep = ir,
                                Method = "Batch-mOpCov",
                                StepSize = NA,
                                N = N,
                                nBatch = Nbatch,
                                Time = difftime(batch_end,batch_start,units = 'secs'),
                                # AV = fit0$av,
                                # EWMABV = fit0$ewmabv,
                                RMSEphi1 = rmseBatch[1],
                                RMSEphi2 = rmseBatch[2],
                                RMSEphi3 = rmseBatch[3],
                                RMSEphi1.avg = rmseBatch[1],
                                RMSEphi2.avg = rmseBatch[2],
                                RMSEphi3.avg = rmseBatch[3]))
      }
      
    }
    
    write.csv(res, filename, row.names = FALSE)
  }
}



