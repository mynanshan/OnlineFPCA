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

idx.start <- 11
idx.end <- 151

sgdtype <- "sgd"



# Plot CIs ---------------------------------------------------------------

seed <- 11
n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("ci_", exprmt, "_sd", seedtext, ".RData")

load(file.path(dirpath, filename))
CIs <- CIobj$CIs
PhiStar <- CIobj$PhiStar
PhiTrue <- CIobj$PhiTrue

idx <- c(51, 101, 151)

pdf(file = file.path(dirpath, "ci1d-ci.pdf"), width = 7, height = 6.2)
pa <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(2, 2, 1, 0))
yranges <- sapply( 1:npc, \(k) range(CIs[,k,,idx]))
for (i in seq_along(idx)) {
  for (k in seq_len(npc)) {
    phik <- PhiTrue[, k]
    phikstar <- PhiStar[, k]
    phik.l <- CIs[, k, 1, idx[i]]
    phik.est <- CIs[, k, 2, idx[i]]
    phik.u <- CIs[, k, 3, idx[i]]
    plot(
      tgrid,
      phikstar,
      type = "l",
      lty = 2,
      lwd = 2,
      xlab = "",
      ylab = "",
      main = NULL,
      ylim = yranges[, k]
    )
    # lines(evalGrid, phik.l)
    # lines(evalGrid, phik.u)
    polygon(
      c(tgrid, rev(tgrid)),
      c(phik.l, rev(phik.u)),
      col = rgb(0.2, 0.3, 0.9, alpha = 0.2),
      border = NA
    )
    # lines(tgrid, phik, lwd = 1.2, lty = 2)
    lines(tgrid, phik.est, lwd = 2, col = 2)
    # lines(evalGrid, phik.l, lty=4, col="#1A73A0", lwd=3)
    if (k == 1) {
      mtext(paste0("n=", (idx[i] - 1) * 100), side = 2, line = 3)
    }
    if (i == 1) {
      mtext(paste0("FPC ", k), side = 3, line = 1, font = 2)
    }
    if (i == length(idx)) {
      mtext("t", side = 1, line = 3, font = 2)
    }
  }
}
par(pa)
dev.off()

# Check coverage rates ---------------------------------------------------

# Current, the coverage rates still deviate observable from the nominal values


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

