library(fda)
library(Matrix)

source("./R/fdaMdim.R")
source("./R/helper.R")

exprmt <- "ci1d"
dirpath <- file.path("experiments", exprmt)

npc <- 3
nrep <- 1000
nt <- 51
t0 <- 0; t1 <- 1
tgrid <- seq(t0,t1,length.out=nt)

nbasis <- 7
basis <- fda::create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = 4)
B <- eval_basis(tgrid, basis)

idx.start <- 11
idx.end <- 151

sgdtype <- "sgd"

cr <- array(0, dim = c(nt, npc, idx.end))

for (seed in seq_len(nrep)) {

  n_digit_seed = ceiling(log10(seed + 1))
  seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
  filename <- paste0("ci_", exprmt, "_sd", seedtext, ".RData")

  errflag <- FALSE
  tryCatch(
    load(file.path(dirpath, filename)),
    error = function(e) {
      cat("Combination unfound:", sgdtype, seed, "\n")
      errflag <<- TRUE
    }
  )
  CIs <- saveObj$CIs[[sgdtype]]

  # optimum
  PhiStar <- as.matrix(B %*% saveObj$sol$Theta)

  # coverage
  for (k in seq_len(npc)) {
    phik_star <- PhiStar[,k]
    inprds <- colMeans(CIs[,k,2,] * phik_star)  # (nt, niter)
    neg.ind <- inprds < 0
    CIs[,k,,neg.ind] <- -CIs[,k,,neg.ind]
    phi_l0 <- CIs[,k,1,]
    phi_u0 <- CIs[,k,3,]
    phi_l <- pmin(phi_l0, phi_u0)
    phi_u <- pmax(phi_l0, phi_u0)
    cr[,k,] <- cr[,k,] + ((phik_star > phi_l) & (phik_star < phi_u)) / nrep
  }
}





read.CI <- function(seed, sgdtype) {

  n_digit_seed = ceiling(log10(seed + 1))
  seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
  filename <- paste0("ci_", exprmt, "_sd", seedtext, ".RData")

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
