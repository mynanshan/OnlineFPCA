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
nPass <- 1
nBlock <- 100
Ninit <- 100
initMethod <- "face"
N <- 5000
stepsize <- 2e-1
stepsize.min <- 5e-2
sgd.step.scale <- 1 # use a smaller step size for sgd
asgd.use <- TRUE
ewmabv.beta.list <- c(0.1, 0.3, 0.5, 0.7, 0.9)

exprmt <- "abv1d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("simu_", exprmt, "_sd", seedtext, ".csv")

set.seed(1234)

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
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)

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

inits <- list(Theta = ThetaInit, lambda = lambdaInit, sigma2 = sigma2Init)

optMethod <- "Adam"

nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)

ewmabv.mat <- c()
for (ewmabv.beta in ewmabv.beta.list) {

  suppressWarnings(
    {
      fit <- fpca.sgd(
        fdata_generator,
        tgrid,
        inits = inits,
        meanfun = FALSE,
        tau = NULL,
        tau.control = list(
          ntau = 1,
          nselect = 1,
          maxtau = 1e-6,
          mintau = 1e-6
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
        ewmabv.beta = ewmabv.beta,
        verbose = FALSE
      )
    }
  )

  vcrits <- fit$vcrit.history

  bv <- vcrits$bv
  ewmabv <- vcrits$ewmabv

  bv <- as.vector(na.omit(c(bv)))
  ewmabv <- as.vector(na.omit(c(ewmabv)))

  ewmabv.mat <- cbind(ewmabv.mat, ewmabv)

}

colnames(ewmabv.mat) <- paste0("beta", ewmabv.beta.list*10)

pdf(file.path(dirpath, paste0("abv_paths.pdf")), width = 7, height = 4.2)
# 1. save current par
oldpar <- par(no.readonly = TRUE)
# 2. bump up right margin so there's room for the legend
#    mar = c(bottom, left, top, right)
par(mar = c(4, 4, 2, 2) + 0.1,
    font.lab = 2)    # bold axis labels
# 3. do your matplot
matplot(
  ewmabv.mat,
  lwd = 2,
  lty = rev(seq_along(ewmabv.beta.list)),
  col = rev(seq_along(ewmabv.beta.list)),
  type = "l",
  xlab = "Number of Data Blocks",
  ylab = "ABV")
# 4. build true plotmath expressions for omega == value
#    Option A: parse() on a character vector
legend_labels <- parse(text = paste0("omega==", ewmabv.beta.list))
# 5. add the legend just **outside** the right side of the plot
legend("bottomright",
       inset   = c(0.05, 0.1),  # negative x‑inset pushes it outside
       xpd     = TRUE,         # allow drawing in the figure margin
       legend  = legend_labels,
       lty     = rev(seq_along(ewmabv.beta.list)),
       col     = rev(seq_along(ewmabv.beta.list)))          # optional: no box around legend
# 6. restore original par
par(oldpar)
dev.off()