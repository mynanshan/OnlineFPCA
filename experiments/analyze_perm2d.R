library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

exprmt <- "perm2d"
dirpath <- file.path("experiments", exprmt)


# FPC UQ Band Plot -------------------------------------------------------

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./data_generation/generator.R")

rp <- get_rp.wang2020(alpha = 2)
m_min <- NULL
m_max <- NULL
m_mean <- 25
m_sd <- 6
m_type <- "gaussian"

t0 <- rp$t0
t1 <- rp$t1
D <- rp$ndim

q <- 3

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

nbasis1d <- 7
basis <- TensorBasis(list(
  create.bspline.basis(c(t0[1], t1[1]), nbasis = nbasis1d, norder = 4),
  create.bspline.basis(c(t0[2], t1[2]), nbasis = nbasis1d, norder = 4)
))
p <- attr(basis, "nbasis")

file_list = list.files(dirpath)
fpc_files = file_list[startsWith(file_list, "Theta")]
nrep = length(fpc_files)

iRecord <- c(10, 30, 50, 100, 150)
nrec <- length(iRecord)
ndata.1rec <- 100
ThetaSGD <- ThetaAdam <- rep(list(array(0, c(p, q, nrep))), nrec)
for (irep in seq_along(fpc_files)) {
  ThetaRecord <- readRDS(file.path(dirpath, fpc_files[irep]))
  for (ip in seq_len(nrec)) {
    ThetaSGD[[ip]][,, irep] <- ThetaRecord[[ip]][['SGD']]
    ThetaAdam[[ip]][,, irep] <- ThetaRecord[[ip]][['Adam']]
  }
}

PhiSGD <- lapply(ThetaSGD, \(Theta) {
  apply(
    Theta,
    3,
    \(X) eval_fd(evalGrid, FuncData(X, basis))
  ) |> # (nt*q, nrep)
    array(c(neval, q, nrep))
})
PhiAdam <- lapply(ThetaAdam, \(Theta) {
  apply(
    Theta,
    3,
    \(X) eval_fd(evalGrid, FuncData(X, basis))
  ) |> # (nt*q, nrep)
    array(c(neval, q, nrep))
})

Phi.CI <- lapply(
  PhiAdam,
  \(Phi) {
    apply(Phi, c(1, 2), quantile, probs = c(0.025, 0.975)) |>
      array(dim = c(2, neval, q)) |>
      aperm(c(2, 3, 1))
  }
)


# install.packages("plot3D")
library(plot3D)

# --- which entries to show and compute their sample sizes
record_indices <- c(1, 3, 5)
sample_sizes <- iRecord[record_indices] * ndata.1rec

# set up a PDF (or PNG) device
pdf(file = file.path(dirpath, "perm2d-ci.pdf"), width = 7, height = 6.2)

# 3×3 panels, with some outer margin for axis labels
par(
  mfrow = c(3, 3),
  mar = c(1, 2, 2, 0),
  oma = c(0, 3, 2, 0),
  xpd = TRUE
)

for (ri in seq_along(record_indices)) {
  idx <- record_indices[ri]

  for (k in seq_len(q)) {
    # reshape back to matrix
    Z_true <- matrix(PhiTrueEval[, k], nevalList[1], nevalList[2])
    Z_lower <- matrix(Phi.CI[[idx]][, k, 1], nevalList[1], nevalList[2])
    Z_upper <- matrix(Phi.CI[[idx]][, k, 2], nevalList[1], nevalList[2])

    # draw the true FPC surface
    persp3D(
      x = matrix(evalGrid[, 1], nevalList[1], nevalList[2]),
      y = matrix(evalGrid[, 2], nevalList[1], nevalList[2]),
      z = Z_true,
      col = rgb(1., 0.8, 0.3, alpha = 0.2),
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
    phi_range <- c(min(Z_true), mean(range(Z_true)), max(Z_true))
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
    surf3D(
      x = matrix(evalGrid[, 1], nevalList[1], nevalList[2]),
      y = matrix(evalGrid[, 2], nevalList[1], nevalList[2]),
      z = Z_lower,
      col = rgb(0.2, 0.3, 0.9, alpha = 0.2),
      border = NA,
      add = TRUE
    )

    # overlay the upper‐CI surface
    surf3D(
      x = matrix(evalGrid[, 1], nevalList[1], nevalList[2]),
      y = matrix(evalGrid[, 2], nevalList[1], nevalList[2]),
      z = Z_upper,
      col = "blue",
      alpha = 0.2,
      border = NA,
      add = TRUE
    )

    # annotate
    if (ri == 1) {
      mtext(paste("FPC", k), side = 3, line = 1, font = 2)
    }
    if (k == 1) {
      mtext(
        paste0("n = ", sample_sizes[ri]),
        side = 2,
        line = 1,
        font = 2
      )
    }
  }
}

dev.off()
