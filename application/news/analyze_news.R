#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fda)
})

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/onlineFPCA.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ================================================================
## User settings
## ================================================================
project_root <- "."
processed_file <- file.path(project_root, "data", "news", "processed", "news_preprocessed.rds")
result_file <- file.path(project_root, "application", "result_news.Rdata")
result_uq_file <- file.path(project_root, "application", "result_news_uq.Rdata")
result_olcov_file <- file.path(project_root, "application", "result_news_olcov.Rdata")
out_dir <- file.path(project_root, "application")


## ================================================================
## Helpers
## ================================================================
clock_label <- function(x) sprintf("%02d:00", as.integer(round(x)) %% 24)

nearest_indices <- function(target, available) {
  vapply(target, function(x) which.min(abs(available - x)), integer(1))
}

get_array_slice <- function(arr, idx) {
  if (length(dim(arr)) != 3L) stop("Expected a 3D array.")
  arr[, , idx, drop = FALSE][, , 1]
}

safe_range <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) c(-1, 1) else range(vals)
}

align_curve_sign <- function(target, reference, basis) {
  target_fd <- smooth.basis(seq(0, 1, length.out = length(target)), target, basis)$fd
  reference_fd <- smooth.basis(seq(0, 1, length.out = length(reference)), reference, basis)$fd
  if (inprod(target_fd, reference_fd) < 0) -target else target
}

align_ci_band <- function(lower, middle, upper, reference, basis) {
  middle_fd <- smooth.basis(seq(0, 1, length.out = length(middle)), middle, basis)$fd
  reference_fd <- smooth.basis(seq(0, 1, length.out = length(reference)), reference, basis)$fd
  if (inprod(middle_fd, reference_fd) < 0) {
    list(lower = -upper, middle = -middle, upper = -lower)
  } else {
    list(lower = lower, middle = middle, upper = upper)
  }
}

match_record_grid <- function(PhiAvgAll, record_subjects, N, nBlock) {
  n_slices <- dim(PhiAvgAll)[3]
  if (length(record_subjects) == n_slices) {
    record_subjects
  } else if (length(record_subjects) == n_slices + 1L) {
    record_subjects[-1]
  } else {
    pmin(seq_len(n_slices) * nBlock, N)
  }
}

## ================================================================
## Load data and fits
## ================================================================
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
news_obj <- readRDS(processed_file)
load(result_file)

has_uq <- file.exists(result_uq_file)
has_compare <- file.exists(result_olcov_file)

if (has_uq) {
  load(result_uq_file)
}
if (has_compare) {
  load(result_olcov_file)
}

par.old <- par(no.readonly = TRUE)

## Prefer the UQ fit for final summaries if available
PhiMainEval <- if (has_uq) PhiAvgEval.uq else PhiAvgEval
PhiMainAll <- if (has_uq) PhiAvgAll.uq else PhiAvgAll
lambda.main <- if (has_uq) lambda.avg.uq else lambda.avg
record_main <- if (has_uq) record_subjects.uq else record_subjects
sgd_time_main <- if (has_uq) sgd_time.uq else sgd_time

evalGrid <- run_config$evalGrid
clock_eval <- run_config$clock_eval
nBlock <- run_config$nBlock
N <- run_config$N
q <- run_config$q
basis <- create.bspline.basis(c(0, 1), nbasis = run_config$nbasis, norder = run_config$basis_order)

record_grid <- match_record_grid(PhiMainAll, record_main, N, nBlock)

## Order FPCs by final eigenvalues
ord <- order(lambda.main, decreasing = TRUE)
PhiMainEval <- PhiMainEval[, ord, drop = FALSE]
PhiMainAll <- PhiMainAll[, ord, , drop = FALSE]
lambda.main <- lambda.main[ord]

if (has_compare) {
  ord.ll <- order(lambdaEst.ll, decreasing = TRUE)
  PhiEst.ll <- PhiEst.ll[, ord.ll, drop = FALSE]
  PhiEstAll.ll <- PhiEstAll.ll[, ord.ll, , drop = FALSE]
  lambdaEst.ll <- lambdaEst.ll[ord.ll]

  q <- min(q, ncol(PhiEst.ll))
  PhiMainEval <- PhiMainEval[, seq_len(q), drop = FALSE]
  PhiMainAll <- PhiMainAll[, seq_len(q), , drop = FALSE]
  lambda.main <- lambda.main[seq_len(q)]
  PhiEst.ll <- PhiEst.ll[, seq_len(q), drop = FALSE]
  PhiEstAll.ll <- PhiEstAll.ll[, seq_len(q), , drop = FALSE]
  lambdaEst.ll <- lambdaEst.ll[seq_len(q)]
  clock_eval_short <- seq(clock_eval[1], tail(clock_eval, 1), length.out = dim(PhiEstAll.ll)[1])
}
if (has_uq) {
  resCI <- resCI[, ord, , , drop = FALSE]
  resCI <- resCI[, seq_len(q), , , drop = FALSE]
}

## Flip signs for readability using the final OnlineFPCA estimate
flip_idx <- colMeans(PhiMainEval) < 0
if (any(flip_idx)) {
  PhiMainEval[, flip_idx] <- -PhiMainEval[, flip_idx, drop = FALSE]
  PhiMainAll[, flip_idx, ] <- -PhiMainAll[, flip_idx, , drop = FALSE]
  if (has_uq) {
    resCI[, flip_idx, , ] <- -resCI[, flip_idx, c(3, 2, 1), , drop = FALSE]
  }
}

fve <- lambda.main / sum(lambda.main)
summary_tbl <- data.frame(
  FPC = seq_len(q),
  eigenvalue = lambda.main[seq_len(q)],
  FVE = fve[seq_len(q)],
  cumFVE = cumsum(fve)[seq_len(q)]
)
write.csv(summary_tbl, file.path(out_dir, "news_fve.csv"), row.names = FALSE)
saveRDS(PhiMainEval, file.path(out_dir, "eigf_news.rds"))
saveRDS(lambda.main, file.path(out_dir, "eigval_news.rds"))

## ================================================================
## Mean plot from preprocessing
## ================================================================
mean_df <- news_obj$mean_curve
pdf(file = file.path(out_dir, "news_mean.pdf"), width = 7, height = 4.8)
plot(
  mean_df$clock, mean_df$mean_raw,
  type = "l", lwd = 3,
  xlab = "Clock time", ylab = "Estimated mean log-click count",
  xaxt = "n", main = NULL
)
axis(1, at = seq(6, 30, by = 4), labels = clock_label(seq(6, 30, by = 4)))
box()
dev.off()

## ================================================================
## Eigenfunction figure
## ================================================================

# Subjects at which to show checkpoint figures
checkpoint_subjects <- c(2000, 5000, 8000)
cp_idx <- nearest_indices(checkpoint_subjects, record_grid)
cp_subjects <- record_grid[cp_idx]
cp_labels <- format(cp_subjects, big.mark = ",")
add_ci <- has_uq
npanel <- length(cp_idx) * q
xlab_panel <- npanel + 1L
legend_panel <- npanel + 2L
ci_sub_id <- 6 * seq(0, length(clock_eval)-1) + 1
pdf(file = file.path(out_dir, "news_eigfun.pdf"), width = 8.2, height = 6.8)
layout(
  mat = matrix(c(seq_len(npanel), rep(xlab_panel, q), rep(legend_panel, q)), byrow = TRUE,
    nrow = length(cp_idx) + 2L, ncol = q
  ),
  heights = c(rep(0.31, length(cp_idx)), 0.04, 0.10)
)
par(mar = c(2.3, 4.3, 2.0, 1.2), oma = c(2, 0, 1, 0))
for (ii in seq_along(cp_idx)) {
  kk <- cp_idx[ii]
  for (jj in seq_len(q)) {
    phik.est <- get_array_slice(PhiMainAll, kk)[, jj]
    phik.compare <- NULL
    if (has_compare) {
      kk2 <- nearest_indices(cp_subjects[ii], record_subjects.ll)[1]
      phik.compare <- get_array_slice(PhiEstAll.ll, kk2)[, jj]
      phik.compare <- align_curve_sign(phik.compare, phik.est, basis)
    }
    yrng <- safe_range(phik.est, phik.compare)
    if (add_ci) {
      ci_idx <- nearest_indices(cp_subjects[ii], match_record_grid(PhiAvgAll.uq, record_subjects.uq, N, nBlock))[1]
      phik.l <- resCI[ci_sub_id, jj, 1, ci_idx]
      phik.m <- resCI[ci_sub_id, jj, 2, ci_idx]
      phik.u <- resCI[ci_sub_id, jj, 3, ci_idx]
      ci_band <- align_ci_band(phik.l, phik.m, phik.u, phik.est, basis)
      phik.l <- ci_band$lower
      phik.m <- ci_band$middle
      phik.u <- ci_band$upper
      yrng <- safe_range(yrng, phik.l, phik.m, phik.u)
    }
    pad <- 0.08 * diff(yrng)
    if (!is.finite(pad) || pad == 0) pad <- 0.3
    yrng <- c(yrng[1] - pad, yrng[2] + pad)
    plot(
      clock_eval, phik.est,
      type = "l", lwd = 3, col = "#A21B2D",
      xlab = "", ylab = "", xaxt = "n",
      ylim = yrng
    )
    axis_labels <- if (ii == length(cp_idx)) clock_label(seq(6, 30, by = 6)) else FALSE
    axis(1, at = seq(6, 30, by = 6), labels = axis_labels)
    if (has_compare) {
      lines(clock_eval_short, phik.compare, lwd = 3, lty = 4, col = "#1A73A0")
    }
    if (add_ci) {
      polygon(
        c(clock_eval, rev(clock_eval)),
        c(phik.l, rev(phik.u)),
        col = rgb(0.2, 0.3, 0.9, alpha = 0.18),
        border = NA
      )
      lines(clock_eval, phik.est, lwd = 3, col = "#A21B2D")
      if (has_compare) {
        lines(clock_eval_short, phik.compare, lwd = 3, lty = 4, col = "#1A73A0")
      }
    }
    if (jj == 1) {
      mtext(paste0("n=", cp_labels[ii]), side = 2, line = 2.8)
    }
    if (ii == 1) {
      mtext(sprintf("FPC %d, %.1f%%", jj, 100 * fve[jj]), side = 3, line = 0.7)
    }
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.53, 0.44, "Clock time", cex = 1.3)
par(mar = c(0, 2, 0.2, 0))
plot.new()
leg_txt <- c("OnlineFPCA-RSGD")
leg_lty <- c(1)
leg_lwd <- c(2.5)
leg_pch <- c(NA)
leg_ptcex <- c(NA)
leg_col <- c("#A21B2D")
if (has_compare) {
  leg_txt <- c(leg_txt, "OnlineCov")
  leg_lty <- c(leg_lty, 4)
  leg_lwd <- c(leg_lwd, 2.5)
  leg_pch <- c(leg_pch, NA)
  leg_ptcex <- c(leg_ptcex, NA)
  leg_col <- c(leg_col, "#1A73A0")
}
# if (add_ci) {
#   leg_txt <- c(leg_txt, "Pointwise CI")
#   leg_lty <- c(leg_lty, NA)
#   leg_lwd <- c(leg_lwd, NA)
#   leg_pch <- c(leg_pch, 15)
#   leg_ptcex <- c(leg_ptcex, 2)
#   leg_col <- c(leg_col, rgb(0.2, 0.3, 0.9, alpha = 0.25))
# }
legend(
  "center",
  legend = leg_txt,
  lty = leg_lty,
  lwd = leg_lwd,
  pch = leg_pch,
  pt.cex = leg_ptcex,
  col = leg_col,
  text.width = max(strwidth(leg_txt)) * 1.18,
  horiz = TRUE, bty = "n", cex = 1.05
)
par(par.old)
dev.off()

## ================================================================
## Tau path
## ================================================================
pdf(file = file.path(out_dir, "news_taupath.pdf"), width = 7.2, height = 4.5)
if (exists("tau.select.uq")) {
  with(tau.select.uq, plot.tau_path2(tau.history, tau.selectId))
} else if (exists("tau.select")) {
  with(tau.select, plot.tau_path2(tau.history, tau.selectId))
} else {
  plot.new()
  text(0.5, 0.5, "Tau-path history not found in the saved result.")
}
dev.off()

## ================================================================
## Time comparison
## ================================================================
pdf(file = file.path(out_dir, "news_time.pdf"), width = 6.5, height = 4.8)
if (has_compare) {
  x1 <- if (length(sgd_time_main) == length(record_grid) + 1L) c(0, record_grid) else record_grid
  y1 <- cumsum(as.numeric(sgd_time_main)) / 60
  matplot(
    x = cbind(x1 / 1000, record_subjects.ll / 1000),
    y = cbind(y1, cumsum(as.numeric(loclin_time)) / 60),
    type = "l", lwd = 3,
    col = c("#A21B2D", "#1A73A0"), lty = c(1, 2),
    xlab = "Number of news items (thousand)",
    ylab = "CPU time (min)"
  )
  legend(
    "topleft",
    legend = c("OnlineFPCA-RSGD", "OnlineCov"),
    lty = c(1, 2), lwd = 2.5,
    col = c("#A21B2D", "#1A73A0"), bty = "n"
  )
} else {
  x1 <- if (length(sgd_time_main) == length(record_grid) + 1L) c(0, record_grid) else record_grid
  plot(
    x1 / 1000, cumsum(as.numeric(sgd_time_main)) / 60,
    type = "l", lwd = 3, col = "#A21B2D",
    xlab = "Number of news items (thousand)",
    ylab = "CPU time (min)"
  )
}
dev.off()

## ================================================================
## A short text summary for the draft
## ================================================================
summary_txt <- file.path(out_dir, "news_summary.txt")
con <- file(summary_txt, open = "wt")
writeLines(c(
  sprintf("Retained morning-news subjects after preprocessing: %d", length(news_obj$dat$Ly)),
  sprintf("Subjects used in the OnlineFPCA run: %d", N),
  sprintf("Dense retained measurements per subject: %d", unique(news_obj$dat$Lmi)[1]),
  sprintf("Checkpoint subjects shown in the eigenfunction figure: %s", paste(cp_labels, collapse = ", ")),
  "",
  "Final OnlineFPCA eigenvalue / FVE summary:",
  paste(capture.output(print(summary_tbl, row.names = FALSE)), collapse = "\n")
), con = con)
close(con)
