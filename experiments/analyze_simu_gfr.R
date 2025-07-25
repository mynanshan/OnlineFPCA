library(dplyr)
library(tidyr)
library(readr)

exprmt <- "gfr"
dirpath <- file.path("experiments", exprmt)


# RMSE plot and table ----------------------------------------------------

res <- do.call(
  rbind,
  lapply(
    1:100, \(i) {
      seed <- i
      n_digit_seed = ceiling(log10(seed + 1))
      seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
      filename <- paste0("simu_",exprmt,"_sd",seedtext,".csv")
      read_csv(file.path(dirpath, filename), show_col_types = FALSE)
    }
  )
)

res <- res |> 
  mutate(
    seed = as.integer(seed),
    Ninit = as.integer(Ninit),
    npc = as.integer(npc),
    nBatch = as.integer(nBatch)
  )

sumres <- res %>% group_by(Method, Ninit, npc, initMethod, N) %>%
  summarise_at(vars(Time, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg), list(mean=mean, sd=sd))

sumres[['SumTime']] <- NA
for (m in unique(sumres$Method)) {
  times <- sumres$Time_mean[sumres$Method==m]
  sumres$SumTime[sumres$Method==m] <- cumsum(times)
}

q <- 3
N0 <- c(5000, 10000, 15000)
N1 = 5000

# zoom-in version
# setEPS()
# postscript(file="text-onlineFPCA/figures/fpca1d_rmse-time_n5000epoch3alpha2.eps",
#            width = 9, height=4)
pa <- par(no.readonly = TRUE)
m <- matrix(c(1,2,3,4,4,4), nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(0.85,0.15))
par(mar=c(5,4,2,2))
for (k in 1:q) {
  # plot onlineFPCA-Adam
  chart.pspladam <- list(
    x = sumres |> filter(Method=="Pspline-Adam") |> pull("SumTime"),
    y = sumres |> filter(Method=="Pspline-Adam") |> pull(paste0("RMSEphi",k,".avg_mean")),
    xm = sumres |> filter(Method=="Pspline-Adam") |> filter(N %in% N0) |>
      arrange(N) |> pull("SumTime"),
    ym = sumres |> filter(Method=="Pspline-Adam") |> filter(N %in% N0) |>
      arrange(N) |> pull(paste0("RMSEphi",k,".avg_mean"))
  )
  chart.psplsgd <- list(
    x = sumres |> filter(Method=="Pspline-SGD") |> pull("SumTime"),
    y = sumres |> filter(Method=="Pspline-SGD") |> pull(paste0("RMSEphi",k,".avg_mean")),
    xm = sumres |> filter(Method=="Pspline-SGD") |> filter(N %in% N0) |>
      arrange(N) |> pull("SumTime"),
    ym = sumres |> filter(Method=="Pspline-SGD") |> filter(N %in% N0) |>
      arrange(N) |> pull(paste0("RMSEphi",k,".avg_mean"))
  )
  chart.loclin <- list(
    x = sumres |> filter(Method=="LocLin") |> pull("SumTime"),
    y = sumres |> filter(Method=="LocLin") |> pull(paste0("RMSEphi",k,".avg_mean")),
    xm = sumres |> filter(Method=="LocLin") |> filter(N == N0[1]) |>
      pull("SumTime"),
    ym = sumres |> filter(Method=="LocLin") |> filter(N == N0[1]) |>
      pull(paste0("RMSEphi",k,".avg_mean"))
  )
  # , xlim=c(0,100)
  pspl.starttime = sumres |>
    filter(Method=="Pspline-SGD") |> pull("SumTime") |> min()
  pspl.endtime = sumres |>
    filter(Method=="Pspline-SGD") |> pull("SumTime") |> max(na.rm = TRUE)
  loclin.starttime = sumres |>
    filter(Method=="LocLin") |> pull("SumTime") |> min()
  loclin.endtime = sumres |>
    filter(Method=="LocLin") |> pull("SumTime") |> max()
  # ymax = sumres |>
  #   pull(paste0("RMSEphi3.avg_mean")) |> max()
  if (k==1) {
    numtext = "1st"
  } else if (k==2) {
    numtext = "2nd"
  } else {
    numtext = "3rd"
  }
  plot(NULL,
       xlim=c(0, pspl.endtime),
       ylim=c(0, 0.6) * c(0.8,1.05),
       xlab="Time (s)", ylab="RMSE", main=paste(numtext, "FPC"))
  # abline(h=chart.loclin$ym, col="#1A73A0", lwd=1, lty=3)
  lines(chart.pspladam$x, chart.pspladam$y, col="#A21B2D", lwd=2, lty=1)
  points(chart.pspladam$xm, chart.pspladam$ym, col="#A21B2D", lwd=1, cex=2, pch=1)
  lines(chart.psplsgd$x, chart.psplsgd$y, col="#0F7323", lwd=2, lty=2)
  points(chart.psplsgd$xm, chart.psplsgd$ym, col="#0F7323", lwd=1, cex=2, pch=2)
  # lines(chart.loclin$x, chart.loclin$y, col="#1A73A0", lwd=2, lty=4)
  # points(chart.loclin$xm, chart.loclin$ym, col="#1A73A0", lwd=1, cex=2, pch=3)
  # legend("topright", xjust=1, legend = c("onlineFPCA-RAdam", "onlineFPCA-RSGD", "onlineCov"),
  #        lty = 1:3, lwd=1.3, pch=1:3, col=c("#A21B2D", "#0F7323", "#1A73A0"), cex=0.8)
}
par(mar=c(1,1,0,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom", inset=c(0.1, 0.2), 
       legend = c("onlineFPCA-RAdam    ", "onlineFPCA-RSGD    "),
       lty = c(1,2), lwd=1.3, pch=1:2, col=c("#A21B2D", "#0F7323"),
       cex=1.1, horiz=TRUE, bty="n")
par(pa)
# dev.off()
# par(mar=c(1,1,0,1))
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# legend("bottom", inset=c(0.1, 0.2), 
#        legend = c("onlineFPCA-RAdam    ", "onlineFPCA-RSGD    ", "onlineCov"),
#        lty = c(1,2,4), lwd=1.3, pch=1:3, col=c("#A21B2D", "#0F7323", "#1A73A0"),
#        cex=1.1, horiz=TRUE, bty="n")
# par(par.old)
# dev.off()


tabres <- ungroup(sumres) |>
  filter(stringr::str_detect(Method, "Batch") | N %in% N0) |>
  dplyr::select(Method, N, SumTime, RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean)
for (i in 1:nrow(tabres)) {
  cat(paste(c(as.character(tabres[i,1:2]), round(as.numeric(tabres[i,3]),1),
              round(as.numeric(tabres[i,4:6]),3)), collapse=" & "), "\\\\ \n")
}

tabres = tabres |> filter(Method != "Batch-PACE")
# PACE is optimized by C codes. Performance too different


meth_ord = c("Pspline-SGD", "Pspline-Adam", "LocLin",
             "Batch-FACE", "Batch-REML", "Batch-SOAP")
tabres = tabres |> 
  mutate(Method = factor(Method, levels=meth_ord)) |> 
  arrange(Method, N)

meth_list = c("Batch-FACE", "Batch-REML", "Batch-SOAP",
              "LocLin", "Pspline-Adam", "Pspline-SGD")
col_list = c("#900B64", "#D77718", "#857F79", "#1A73A0", "#A21B2D", "#0F7323")
pch_list = c(0, 3, 4, 5, 1, 2)
lab_list = c("FACE", "REML", "SOAP", "onlineCov",
             "onlineFPCA-RAdam", "onlineFPCA-RSGD")
tabres[['Color']] = col_list[match(tabres$Method, meth_list)]
tabres[['Shape']] = pch_list[match(tabres$Method, meth_list)]
tabres[['Label']] = lab_list[match(tabres$Method, meth_list)]
time_range = range(tabres$SumTime)
rmse_range = range(as.matrix(tabres[,paste0('RMSEphi', 1:3, '.avg_mean')]) |> c())


setEPS()
postscript("text-onlineFPCA/figures/fpca1d-sim.eps", width=11.4, height=4.2)
layout(mat = matrix(c(1,2,3,4), nrow=1, byrow=TRUE),
       widths = c(0.28,0.28,0.28,0.16))
par(mar=c(5,3,3,2), oma=c(0,3,0,0))
for (k in 1:q) {
  curr_var = paste0('RMSEphi', k, '.avg_mean')
  # if (k==1) {
  #   par(oma=c(4,0,0,0))
  # } else {
  #   par(oma=c(0,0,0,0))
  # }
  if (k==1) {
    numtext = "1st"
  } else if (k==2) {
    numtext = "2nd"
  } else {
    numtext = "3rd"
  }
  plot(NULL, xlim=c(0, max(time_range)), ylim=rmse_range,
       xlab="CPU Time (s)", ylab="", main=paste(numtext, "FPC"))
  if (k==1) mtext("RMSE", side=2, line=3)
  points(tabres |> pull(SumTime), tabres |> pull(curr_var), cex=2, lwd=2,
        col = tabres |> pull(Color), pch = tabres |> pull(Shape))
}
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("left",
       inset=0,
       legend = unique(tabres$Label), lty=0,
       lwd=2, pch=unique(tabres$Shape),
       col=unique(tabres$Color),
       cex=1.3, horiz=FALSE, bty="n")
par(par.old)
dev.off()



# FPC plots --------------------------------------------------------------

source("./R/fdaMdim.R")
source("./data_generation/generator.R")

dat <- gendata.gfr(1)
tgrid <- evalGrid <- dat$tgrid
t0 <- dat$t0
t1 <- dat$t1
neval <- length(evalGrid)
muTrueEval <- numeric(neval)
PhiTrueEval <- dat$Phi
lambdaTrue <- dat$lam
q <- npc <- length(lambdaTrue)

nbasis <- 7
basis <- create.bspline.basis(c(t0,t1), nbasis = nbasis, norder = 4)
p <- nbasis

file_list = list.files(dirpath)
fpc_files = file_list[startsWith(file_list, "Theta")]
nrep = length(fpc_files)

nPass <- 3
ThetaSGD <- ThetaAdam <- rep(list(matrix(0, nrow = p, ncol = q)), 3)
for (f in fpc_files) {
  ThetaRecord <- readRDS(file.path(dirpath, f))  
  for (ip in seq_len(nPass)) {
    ThetaSGD[[ip]] <- ThetaSGD[[ip]] + ThetaRecord[[ip]][['SGD']] / nrep
    ThetaAdam[[ip]] <- ThetaAdam[[ip]] + ThetaRecord[[ip]][['Adam']] / nrep
  }
}

PhiSGD <- lapply(ThetaSGD, \(Theta) eval_fd(evalGrid, FuncData(Theta, basis)))
PhiAdam <- lapply(ThetaAdam, \(Theta) eval_fd(evalGrid, FuncData(Theta, basis)))


dat.true <- data.frame(
  irep=0,
  Method="Truth",
  t=evalGrid,
  phi1=PhiTrueEval[,1],
  phi2=PhiTrueEval[,2],
  phi3=PhiTrueEval[,3]
)
dat.sgd <- do.call(
  bind_rows,
  lapply(
    seq_along(PhiSGD), \(irep) {
      data.frame(
        irep=irep,
        Method="OnlineFPCA-SGD",
        t=evalGrid,
        phi1=PhiSGD[[nPass]][,1,irep],
        phi2=PhiSGD[[nPass]][,2,irep],
        phi3=PhiSGD[[nPass]][,3,irep]
      )
    }
  )
)
dat.sgd <- data.frame(
  Method="OnlineFPCA-SGD",
  t=evalGrid,
  phi1=PhiSGD[[nPass]][,1],
  phi2=PhiSGD[[nPass]][,2],
  phi3=PhiSGD[[nPass]][,3]
)
dat.adam <- data.frame(
  Method="OnlineFPCA-Adam",
  t=evalGrid,
  phi1=PhiAdam[[nPass]][,1],
  phi2=PhiAdam[[nPass]][,2],
  phi3=PhiAdam[[nPass]][,3]
)
plotdat <- bind_rows(dat.true, dat.sgd, dat.adam) %>%
  pivot_longer(
    cols = starts_with("phi"),
    names_to  = "phiId",
    values_to = "phiVal"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("Truth", "OnlineFPCA-SGD", "OnlineFPCA-Adam")),
    # map “phi1”→“phi[1]” etc so label_parsed knows to render φ₁, φ₂, φ₃
    phiId = recode(phiId,
      phi1 = "phi[1](t)",
      phi2 = "phi[2](t)",
      phi3 = "phi[3](t)"
    )
  )

ggplot(plotdat) +
  geom_line(aes(x = t, y = phiVal)) +
  facet_grid(
    rows     = vars(phiId),
    cols     = vars(Method),
    labeller = label_parsed,     # interpret phi[1], phi[2], phi[3] as expressions
    scales = "free_y"
  ) +
  theme_bw() +
  theme(
    panel.grid    = element_blank(),
    axis.text.x = element_text(size=8)
  ) +
  labs(
    x    = expression(t),    # math axis labels
    y    = "FPC Value",
    fill = NULL                 # drop legend title
  ) +
  scale_fill_viridis_c()        # continuous viridis palette

ggsave(
  file.path(dirpath, "simu-gfr-fpc.pdf"), device = "pdf",
  width = 5.4, height = 4.8
)
