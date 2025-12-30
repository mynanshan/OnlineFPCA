library(fda)
library(Matrix)

source("./R/fdaMdim.R")
source("./R/helper.R")

exprmt <- "ci1d"
dirpath <- file.path("experiments", exprmt)

npc <- 3
nt <- 51
t0 <- 0
t1 <- 1
tgrid <- seq(t0, t1, length.out = nt)

idx.start <- 12
idx.end <- 151

sgdtype <- "sgd"


# Plot CIs ---------------------------------------------------------------

seed <- 2
n_digit_seed <- ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("ci_", exprmt, "_sd", seedtext, ".RData")

load(file.path(dirpath, filename))
CIs <- CIobj$CIs
PhiStar <- CIobj$PhiStar
PhiTrue <- CIobj$PhiTrue

idx <- c(21, 36, 51)

pdf(file = file.path(dirpath, "ci1d-ci.pdf"), width = 7, height = 6.2)
pa <- par(no.readonly = TRUE)
share_y <- FALSE
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(2, 2, 1, 0))
yranges <- sapply(1:npc, \(k) range(CIs[, k, , idx]))
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
      ylim = if (share_y) range(yranges) else yranges[, k]
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
