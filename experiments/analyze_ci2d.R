library(fda)
library(Matrix)

source("./R/fdaMdim.R")
source("./R/helper.R")

exprmt <- "ci2d"
dirpath <- file.path("experiments", exprmt)

npc <- 3
nt <- 51
t0 <- 0; t1 <- 1
tgrid <- seq(t0,t1,length.out=nt)

nevalList <- c(51, 51)
neval <- prod(nevalList)
t0 <- c(0, 0); t1 <- c(1, 1)
evalGridList <- lapply(1:2, \(d) {
  seq(t0[d], t1[d], length.out = nevalList[d])
})
evalGrid <- margins2grid(evalGridList)

sgdtype <- "sgd"

idx.start <- 12
idx.end <- 251


# Plot CIs ---------------------------------------------------------------

library(plot3D)

seed <- 1
n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("ci_", exprmt, "_sd", seedtext, ".RData")

load(file.path(dirpath, filename))
CIs <- CIobj$CIs
PhiStar <- match_fpc(CIobj$PhiStar, PhiTrue)
PhiTrue <- CIobj$PhiTrue

idx <- c(51, 151, 251)

pdf(file = file.path(dirpath, "ci2d-ci.pdf"), width = 7, height = 6.2)
pa <- par(no.readonly = TRUE)
par(
  mfrow = c(3, 3),
  mar = c(1, 2, 2, 0),
  oma = c(0, 3, 2, 0),
  xpd = TRUE
)
subId <- lapply(nevalList, \(N) {
  tmp <- logical(N)
  tmp[seq(1, N, 2)] <- TRUE
  tmp
})
nevalSub <- sapply(subId, sum)
subId2 <- c(outer(subId[[1]], subId[[2]], "&"))
for (i in idx) {
  for (k in seq_len(npc)) {
    phik <- PhiTrue[, k]
    phikstar <- PhiStar[, k]
    phik.l <- CIs[, k, 1, i]
    phik.est <- CIs[, k, 2, i]
    phik.u <- CIs[, k, 3, i]
    # reshape back to matrix
    z <- matrix(phik[subId2], nevalSub[1], nevalSub[2])
    z.est <- matrix(phik.est[subId2], nevalSub[1], nevalSub[2])
    z.star <- matrix(phikstar[subId2], nevalSub[1], nevalSub[2])
    z.l <- matrix(phik.l[subId2], nevalSub[1], nevalSub[2])
    z.u <- matrix(phik.u[subId2], nevalSub[1], nevalSub[2])
    # draw the true FPC surface
    persp3D(
      x = matrix(evalGrid[subId2, 1], nevalSub[1], nevalSub[2]),
      y = matrix(evalGrid[subId2, 2], nevalSub[1], nevalSub[2]),
      z = z.star,
      # z = z,
      col = rgb(1., 0.8, 0.3, alpha = 0.7),
      border = NA,
      xlab = "",
      ylab = "",
      zlab = "",
      theta = 40,
      phi = 30, # view angles
      expand = 0.6,
      lighting = TRUE
    )
    # text3D for axis labels
    t1_range <- c(
      min(evalGridList[[1]]),
      mean(range(evalGridList[[1]])),
      max(evalGridList[[1]])
    )
    t1_span <- t1_range[3] - t1_range[1]
    t2_range <- c(
      min(evalGridList[[2]]),
      mean(range(evalGridList[[2]])),
      max(evalGridList[[2]])
    )
    t2_span <- t2_range[3] - t2_range[2]
    phi_range <- c(min(z.star), mean(range(z.star)), max(z.star))
    phi_span <- phi_range[3] - phi_range[1]
    text_shift <- 0.2
    # t₁ on the x-axis
    text3D(
      x = t1_range[2],
      y = t2_range[1] - text_shift * t2_span,
      z = phi_range[1] - text_shift * phi_span,
      expression(t[1]),
      add = TRUE
    )
    # t₂ on the y-axis
    text3D(
      x = t1_range[3] + 0.5 * text_shift * t1_span,
      y = t2_range[2],
      z = phi_range[1] - text_shift * phi_span,
      expression(t[2]),
      add = TRUE
    )
    # φₖ on the z-axis
    text3D(
      x = t1_range[1] - text_shift * t1_span,
      y = t2_range[1] - text_shift * t2_span,
      z = phi_range[2],
      as.expression(bquote(phi[.(k)])),
      add = TRUE
    )
    # overlay the lower‐CI surface
    lbound_col <- ubound_col <- rgb(0.6, 0.6, 0.9, alpha = 0.5)
    surf3D(
      x = matrix(evalGrid[subId2, 1], nevalSub[1], nevalSub[2]),
      y = matrix(evalGrid[subId2, 2], nevalSub[1], nevalSub[2]),
      z = z.l,
      col = lbound_col,
      border = NA,
      add = TRUE
    )
    # overlay the upper‐CI surface
    surf3D(
      x = matrix(evalGrid[subId2, 1], nevalSub[1], nevalSub[2]),
      y = matrix(evalGrid[subId2, 2], nevalSub[1], nevalSub[2]),
      z = z.u,
      col = ubound_col,
      alpha = 0.2,
      border = NA,
      add = TRUE
    )
    # annotate
    if (i == idx[1]) {
      mtext(paste("FPC", k), side = 3, line = 1, font = 2)
    }
    if (k == 1) {
      mtext(
        paste0("n = ", (i - 1) * 100),
        side = 2,
        line = 1,
        font = 2
      )
    }
  }
}
par(pa)
dev.off()
