library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

exprmt <- "boot1d"
dirpath <- file.path("experiments", exprmt)


# FPC UQ Band Plot -------------------------------------------------------

source("./R/fdaMdim.R")
source("./data_generation/generator.R")

alpha <- 2
rp <- get_rp.yang2021(alpha = alpha, npc = 10)
t0 <- rp$t0
t1 <- rp$t1
D <- rp$ndim

q <- npc <- 3

# evaluate true mean and eigenfunctions
neval <- 101
evalGrid <- seq(t0, t1, length.out = neval)
muTrueEval <- rp$meanfun(evalGrid)
PhiTrueEval <- rp$eigfun(evalGrid)
lambdaTrue <- rp$eigval

nbasis <- 7
basis <- create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = 4)
p <- nbasis

file_list = list.files(dirpath)
fpc_files = file_list[startsWith(file_list, "Theta")]
nrep = length(fpc_files)

iRecord <- c(10, 30, 50, 100)
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

get.ci <- function(x, method = c("quantile", "normal")) {
  method = match.arg(method)
  if (method == "quantile") {
    return(quantile(x, probs = c(0.025, 0.975)))
  }
  if (method == "normal") {
    z = qnorm(0.975)
    return(mean(x) + c(-1, 1) * z * sd(x))
  }
}

Phi.CI <- lapply(
  PhiAdam,
  \(Phi) {
    apply(Phi, c(1, 2), get.ci, method="normal") |>
      array(dim = c(2, neval, npc)) |>
      aperm(c(2, 3, 1))
  }
)


pdf(file = file.path(dirpath, "boot1d-ci.pdf"), width = 7, height = 6.2)
pa <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(2, 2, 1, 0))
yranges <- do.call(
  rbind,
  lapply(Phi.CI, \(Phi) rbind(apply(Phi, 2, min), apply(Phi, 2, max)))
) |>
  apply(2, range)
for (irec in seq_along(iRecord)) {
  if (irec == 2) {
    next()
  }
  for (k in seq_len(npc)) {
    phik <- PhiTrueEval[, k]
    phik.l <- Phi.CI[[irec]][, k, 1]
    phik.u <- Phi.CI[[irec]][, k, 2]
    plot(
      evalGrid,
      phik,
      type = "l",
      lwd = 2,
      xlab = "",
      ylab = "",
      main = NULL,
      ylim = yranges[, k]
    )
    # lines(evalGrid, phik.l)
    # lines(evalGrid, phik.u)
    polygon(
      c(evalGrid, rev(evalGrid)),
      c(phik.l, rev(phik.u)),
      col = rgb(0.2, 0.3, 0.9, alpha = 0.2),
      border = NA
    )
    lines(evalGrid, phik, lwd = 2)
    # lines(evalGrid, phik.l, lty=4, col="#1A73A0", lwd=3)
    if (k == 1) {
      mtext(paste0("n=", iRecord[irec] * ndata.1rec), side = 2, line = 3)
    }
    if (irec == 1) {
      mtext(paste0("FPC ", k), side = 3, line = 1, font = 2)
    }
    if (irec == length(iRecord)) {
      mtext("t", side = 1, line = 3, font = 2)
    }
  }
}
par(pa)
dev.off()
