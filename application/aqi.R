library(tidyr)
library(dplyr)
library(readr)
library(stringr)
# library(ggplot2)
library(lubridate)
# library(giscoR)
# library(sf)
# library(spData)

library(fda)
library(Matrix)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./R/onlineFPCA.R")
source("./R/fpcaReg.R")

par.old <- par(no.readonly = TRUE)


datapath <- "data/epa-aqs"
respath <- "application"

lonrange <- c(-125, -67)
latrange <- c(25, 49)
yearrange <- 1982:2022

load(file.path(datapath, "aqi-us.Rda"))
siteinfo <- readxl::read_excel(file.path(datapath, "siteinfo.xlsx"))


# dat <- dat |>
#   left_join(siteinfo |>
#               dplyr::select(Latitude, Longitude, Latitude.Binned, Longitude.Binned),
#             by = c("Latitude", "Longitude"))

dat = dat |> drop_na()
sites <- dat |> distinct(Latitude.Binned, Longitude.Binned)

# rescaling latitude and longitude
sites$s <- (sites$Latitude.Binned - latrange[1]) / diff(latrange)
sites$t <- (sites$Longitude.Binned - lonrange[1]) / diff(lonrange)
dat$s <- (dat$Latitude.Binned - latrange[1]) / diff(latrange)
dat$t <- (dat$Longitude.Binned - lonrange[1]) / diff(lonrange)

# log transform of AQI
dat$y <- log1p(dat$AQI)
# par(mfrow=c(1,2))
# hist(dat$AQI, xlab="AQI", main="", breaks = 30)
# hist(dat$y, xlab="log(AQI+1)", main="", breaks = 30)
# par(mfrow=c(1,1))

# Settings ==============================
latgrid <- seq(latrange[1] + 1, latrange[2], 2)
longrid <- seq(lonrange[1] + 1, lonrange[2], 2)
locGrid <- margins2grid(list(latgrid, longrid))
# save(latrange, lonrange, yearrange, locGrid, evalGrid, basis, file="aqi_settings.Rdata")
colnames(locGrid) <- c("lat", "lon")

locGridRescale = scale(
  locGrid,
  center = c(latrange[1], lonrange[1]),
  scale = c(diff(latrange), diff(lonrange))
)

evalGridList <- list(
  (latgrid - latrange[1]) / diff(latrange),
  (longrid - lonrange[1]) / diff(lonrange)
)
evalGrid <- margins2grid(evalGridList)

basis <- TensorBasis(list(
  create.bspline.basis(c(0, 1), nbasis = 6, norder = 4),
  create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
))
p <- attr(basis, "nbasis")

B <- eval_basis(locGridRescale, basis)
G <- get_basis_inprod_matrix(basis)
GR <- chol(G)
invG <- chol2inv(GR)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv

# numbering the locations
lat_id = match(dat$Latitude.Binned, latgrid)
lon_id = match(dat$Longitude.Binned, longrid)
dat$LocId = (lon_id - 1) * length(latgrid) + lat_id


# mean estimation ====================
fit.mean <- mgcv::bam(y ~ te(s, t), data = dat)
sites$y.mean <- mu <- as.numeric(predict(fit.mean, newdata = sites))

# ggplot() +
#   geom_sf(data=us_states) +
#   geom_point(aes(x=Longitude.Binned, y=Latitude.Binned, color=y.mean,
#                  size=y.mean), stroke=1,
#              data=sites, shape=1) +
#   theme_bw() +
#   scale_color_continuous(type = "viridis") +
#   labs(x="Longitude", y="Latitude",
#        color = "Average log(AQI+1)", size = "Average log(AQI+1)")

dat$z <- dat$y - as.numeric(predict(fit.mean, newdata = dat))

Lmi <- dat |> group_by(Date.Local) |> summarise(Nobs = n()) |> pull(Nobs)

# algo settings ==================
N = length(unique(dat$Date.Local))
q <- 6
nBatch <- 10
nParams <- 6
nPass <- 10
nBlock <- 100
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
# stepsizeList <- c(1e0, 3e-1, 1e-1, 3e-2)
stepsize0 <- 0.1
stepsize.min <- 1e-5
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
nBlockIter <- nBlock / nBatch
nIter1pass <- round(N / nBatch)
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)


# # Initialization ====================

# # # # read a sampled dataset
# dat$subjId <- rep(1:length(Lmi), Lmi)
#
# sample_ratio <- 0.1
# Lmi.sub <- if_else(Lmi < 10, Lmi, 10)
#
# Ninit <- 100
# set.seed(2501)
# subjId.sub <- sort(sample(1:length(Lmi), Ninit, replace = FALSE))
#
# dat.sub <- dat |> filter(subjId %in% subjId.sub)
# ind <- lapply(subjId.sub, \(i) {
#   ind <- rep(F, Lmi[i])
#   idx <- sample(1:Lmi[i], Lmi.sub[i], replace = FALSE)
#   ind[idx] <- TRUE
#   return(ind)
# }) |> unlist()
# dat.sub <- dat.sub[ind,]
#
# dat2list = function(dat, subjId="Date.Local") {
#   unique_subj = unique(dat[[subjId]])
#   Ly = Ltid = rep(list(NULL), length(unique_subj))
#   Lmi = numeric(length(unique_subj))
#   for (i in seq_along(unique_subj)) {
#     idx = dat[[subjId]]==unique_subj[i]
#     Ly[[i]] = dat$z[idx]
#     Ltid[[i]] = dat$LocId[idx]
#     Lmi[i] = sum(idx)
#   }
#   return(list(Ly=Ly, Ltid=Ltid, Lmi=Lmi))
# }
#
# # # cov est function
# mOpCov_path <- "external_codes/mOpCov/"
# source(paste0(mOpCov_path, "mOpCov_prep.R"))
# Rcpp::sourceCpp(paste0(mOpCov_path, "mOpCov_cpp.cpp"))
#
# locations <- dat.sub |> dplyr::select(s, t) |> as.matrix()
# responses <- dat.sub$z
# subject_ids <- rep(1:Ninit, Lmi.sub[subjId.sub])
#
# start_init = Sys.time()
# CovRes <- mOpCov(location=locations, x=responses, subject=subject_ids,
#                  q=c(6,6), lam = list(lam = 1e-10, alpha = 1e-6), ker = "cos")
# FpcaOut <- fpca.mOpCov(OUT = CovRes)
# end_init = Sys.time()
# lambdaInit <- FpcaOut$Eigen$values[1:q]
# PhiInitEvalFull <- computeEigen(evalGrid,CovRes,FpcaOut)
# PhiInitEval <- PhiInitEvalFull[,1:q]
#
# time_init = difftime(end_init, start_init)
# save(PhiInitEvalFull, PhiInitEval, lambdaInit, time_init, file=file.path("application/aqi_init.Rdata"))
load(file.path("application/aqi_init.Rdata"))
ThetaInit <- smooth_basis(evalGrid, PhiInitEval, basis, lambda = 1e-8)$coefs
norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- eval_fd(evalGrid, FuncData(ThetaInit, basis))
sigma2Init <- 1e-3

# Online Estimation ====================
date_unique = unique(dat$Date.Local)
subj_start_id = c(1, cumsum(Lmi[-N]) + 1)
subj_end_id = cumsum(Lmi)

fdata_generator <- function(n, total_count) {
  idx <- (total_count):(total_count + n - 1) %% N + 1
  curr_Ly = lapply(idx, \(i) dat$z[subj_start_id[i]:subj_end_id[i]])
  curr_Ltid = lapply(idx, \(i) dat$LocId[subj_start_id[i]:subj_end_id[i]])
  return(list(Ly = curr_Ly, Ltid = curr_Ltid, Lmi = Lmi[idx]))
}

ThetaInit <- manifold.Stiefel.retract(ThetaInit, NULL, G)
grad_Theta_init <- objfun(
  fdata_generator(Ninit,0)$Ly[1:Ninit],
  fdata_generator(Ninit,0)$Ltid[1:Ninit],
  ThetaInit,
  lambdaInit,
  sigma2Init,
  NULL,
  0,
  "grad"
)$grad_Theta
grad_Theta_init <- manifold.Stiefel.project(grad_Theta_init, ThetaInit, G)
g_Theta_init_norm <- sqrt(sum(grad_Theta_init * G %*% grad_Theta_init) / q)
sgd_lr0 <- stepsize0 / g_Theta_init_norm

inits <- list(
  Theta = ThetaInit,
  lambda = pmax(lambdaInit, lambdaInit[1] / (1:q)),
  sigma2 = sigma2Init
)


# tmp = dat2list(dat.sub)
# lik0 = objfun(tmp$Ly, tmp$Ltid, ThetaInit, lambdaInit,
#               sigma2=0.1, theta_mu=NULL, tau=0, stats = "loss")$lik
# pen = sum(ThetaInit * Omega %*% ThetaInit)
# taumax = ceiling(log10(lik0 / pen))
taumax = 1e-3

sgdtype <- "adagrad"

fit <- fpca.sgd(
  fdata_generator,
  tgrid,
  inits = inits,
  meanfun = FALSE,
  npc = 4,
  tau = NULL,
  tau.control = list(
    ntau = nParams,
    nselect = 2,
    maxtau = 10^taumax,
    mintau = 10^(taumax - nParams + 1)
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
  nIter.tauNoIncrease = floor(0.7 * nIter1pass),
  period.tune = nBlockIter,
  period.record = nBlockIter,
  verbose = TRUE
)


save(fit, time_init, file = file.path(respath, "fit_aqi.Rdata"))


## === compare to some covariance estimation ==========

# YBatch <- dat$z
# TBatch <- as.matrix(dat |> dplyr::select(s, t))
# # covariance estimation
# tmpdat = df2Ldata(TBatch, YBatch, Lmi)
# covdat <- cov.form_covdat(tmpdat$Lt, tmpdat$Ly)
# bam_start = Sys.time()
# fit.cov <- mgcv::bam(z ~ t2(s1, s2, t1, t2), data = covdat)
# end_start = Sys.time()
# covdat <- cov.form_covdat(evalGridSmall, diag.rm = FALSE)
# cov.hat <- as.numeric(predict(fit.cov, covdat))
# cov.hat <- matrix(cov.hat,sqrt(length(cov.hat)))
# cov.hat <- cov.symmetrize(cov.hat)
# cov.true <- rp$eigfun(evalGridSmall) %*% diag(rp$eigval) %*% t(rp$eigfun(evalGridSmall))
# image(as((cov.true - cov.hat)[seq(1,440,10),seq(1,440,10)], "dMatrix"))
# # NOTE: no guarantee that cov.hat is positive defBatche
# eigObj <- eigen(cov.hat, symmetric = TRUE)
# eigvecs <- eigObj$vectors[,1:q]
# ThetaBatch <- smooth_basis(evalGridSmall, eigvecs, basis, lambda = 1e-10)$coefs
# lambdaBatch <- eigObj$values[1:q]
# norm_factor <- diag(t(ThetaBatch) %*% G %*% ThetaBatch)
# ThetaBatch <- sweep(ThetaBatch, 2, sqrt(norm_factor), "/")
# lambdaBatch <- lambdaBatch * norm_factor
# PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
# PhiBatchEval <- flip_direc(PhiBatchEval, PhiTrueEval[,1:q])
# flipId <- attributes(PhiBatchEval)$flipped
# ThetaBatch[,flipId] <- -ThetaBatch[,flipId]
# rm(covdat, cov.hat, cov.true, eigObj)
# covdat <- data.frame(cbind(TBatch, TBatch), YBatch^2)
# names(covdat) <- c(paste0("s", 1:D), paste0("t", 1:D), "z")
# covdiag.hat <- as.numeric(predict(fit.cov, covdat))
# sigma2Batch <- mean(covdat$z - covdiag.hat)
# rm(fit.cov, covdat, covdiag.hat)
# if (sigma2Batch <= 1e-3) sigma2Batch <- 1e-3
#
# fit.cov = stats::loess(z ~ s1 + s2 + t1 + t2, data = covdat)

## === compare to soap 2d ==========

batch_start <- Sys.time()
fitBatch <- fpca.reg(
  fdata_generator(N, 0)$Ly,
  fdata_generator(N, 0)$Ltid,
  inits = list(
    Theta = ThetaInit,
    lambda = lambdaInit,
    sigma2 = sigma2Init
  ),
  meanfun = FALSE,
  npc = 6,
  maxIter = 300,
  nu = 0.5,
  verbose = TRUE,
  record_iterations = TRUE,
  tau = 10^seq(-6, -6, 1),
  use_validation_set = FALSE,
  refine_alpha = FALSE
)
batch_end <- Sys.time()
batch_time = difftime(batch_end, batch_start) # 35.40681 mins

save(fitBatch, batch_time, file = file.path(respath, "fit_aqi_batch.Rdata"))
