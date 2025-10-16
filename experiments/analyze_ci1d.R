library(fda)
library(Matrix)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./data_generation/generator.R")

N <- 5000

exprmt <- "ci1d"
dirpath <- file.path("experiments", exprmt)

rp <- get_rp.yang2021(alpha = 2, npc = 10)
m_min <- NULL
m_max <- NULL
m_mean <- 6
m_sd <- 2
m_type <- "gaussian"

t0 <- rp$t0
t1 <- rp$t1
D <- rp$ndim

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
nt <- length(tgrid)

muTrueGrid <- rp$meanfun(tgrid)
PhiTrueGrid <- rp$eigfun(tgrid)
lambdaTrue <- rp$eigval

npc <- 3

sgdtype <- "sgd"

idx.start <- 11
idx.end <- 151

read.CI <- function(seed, sgdtype) {

  n_digit_seed = ceiling(log10(seed + 1))
  seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
  filebasename <- paste0("ci_", exprmt, "_sd", seedtext, "_")

  errflag <- FALSE
  tryCatch(
    CIs <- readRDS(file.path(dirpath, paste0(filebasename, sgdtype, ".rds"))),
    error = function(e) {
      cat("Combination unfound:", sgdtype, seed, "\n")
      errflag <<- TRUE
    }
  )
  if (errflag) return(NULL)

  for (j in seq(idx.start, idx.end)) {
    PhiEstFlipped <- match_fpc(CIs[,,2,j], PhiTrueGrid[,1:npc])
    CIs[,,,j] <- CIs[,attributes(PhiEstFlipped)$match_id,,j]
    CIs[,attributes(PhiEstFlipped)$flipped,,j] <-
      -CIs[,attributes(PhiEstFlipped)$flipped,,j]
  }

  CIs <- CIs[,,,idx.start:idx.end]
  CIs
}


CI.to.cvrg <- function(CIs, PhiTrue, npc = 3) {
  nt <- dim(CIs)[1]
  cvrg <- array(dim = c(nt, npc, idx.end - idx.start + 1))
  for (k in seq_len(npc)) {
    phi_true <- PhiTrueGrid[,k]
    inprds <- colMeans(CIs[,k,2,] * PhiTrueGrid[,k])  # (nt, niter)
    neg.ind <- inprds < 0
    CIs[,k,,neg.ind] <- -CIs[,k,,neg.ind]
    phi_l0 <- CIs[,k,1,]
    phi_u0 <- CIs[,k,3,]
    phi_l <- pmin(phi_l0, phi_u0)
    phi_u <- pmax(phi_l0, phi_u0)
    # bias correction:
    cvrg[,k,] <- (phi_true > phi_l) & (phi_true < phi_u)
  }
  cvrg
}

nrep <- 1000

cvrg.rate <- array(0, dim = c(nt, npc, idx.end - idx.start + 1))

for (ir in seq_len(nrep)) {
  seed <- ir
  n_digit_seed = ceiling(log10(seed + 1))
  seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
  filebasename <- paste0("ci_", exprmt, "_sd", seedtext, "_")
  CIs <- read.CI(ir, sgdtype)
  if (is.null(CIs)) next
  cvrg <- CI.to.cvrg(CIs, PhiTrueGrid, npc)
  cvrg.rate <- cvrg.rate + cvrg / nrep
}

CI.to.ocid <- function(CIs, PhiTrue, npc = 3) {
  nt <- dim(CIs)[1]
  ocid <- array(dim = c(nt, npc, idx.end - idx.start + 1))
  for (k in seq_len(npc)) {
    phi_true <- PhiTrueGrid[,k]
    inprds <- colMeans(CIs[,k,2,] * PhiTrueGrid[,k])  # (nt, niter)
    neg.ind <- inprds < 0
    CIs[,k,,neg.ind] <- -CIs[,k,,neg.ind]
    phi_l0 <- CIs[,k,1,]
    phi_u0 <- CIs[,k,3,]
    phi_l <- pmin(phi_l0, phi_u0)
    phi_u <- pmax(phi_l0, phi_u0)
    # bias correction:
    phi_l <- phi_l - bias[,k,]
    phi_u <- phi_u - bias[,k,]
    inCIind <- (phi_true > phi_l) & (phi_true < phi_u)
    ocid[,k,] <- inCIind * 0 +
      (1 - inCIind) * (
        pmin(phi_true - phi_l, 0) + pmax(phi_true - phi_u, 0)
      )
  }
  ocid
}

ocid <- array(0, dim = c(nt, npc, idx.end - idx.start + 1, nrep))

for (ir in seq_len(nrep)) {
  seed <- ir
  n_digit_seed = ceiling(log10(seed + 1))
  seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
  filebasename <- paste0("ci_", exprmt, "_sd", seedtext, "_")
  CIs <- read.CI(ir, sgdtype)
  if (is.null(CIs)) next
  ocid[,,,ir] <- CI.to.ocid(CIs, PhiTrueGrid, npc)
}
CIs <- read.CI(10, sgdtype)


# Check bias -------------------------------------------------------------

PhiEstMean <- array(0, dim = c(nt, npc, idx.end - idx.start + 1))

for (ir in seq_len(nrep)) {
  seed <- ir
  n_digit_seed = ceiling(log10(seed + 1))
  seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
  filebasename <- paste0("ci_", exprmt, "_sd", seedtext, "_")
  CIs <- read.CI(ir, sgdtype)
  if (is.null(CIs)) next

  PhiEstMean <- PhiEstMean + CIs[,,2,] / nrep
}

bias <- sweep(PhiEstMean, c(1,2), PhiTrueGrid[,1:npc], "-")
bias_norm <- apply(
  PhiEstMean, 3, \(Phi) {
    colMeans((Phi - PhiTrueGrid[,1:npc])^2)
  }
)
