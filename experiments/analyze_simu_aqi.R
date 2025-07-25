library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)

exprmt <- "aqi"
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
    nBatch = as.integer(nBatch)
  )

sumres <- res %>% group_by(Method, N) %>%
  summarise_at(vars(
    Time, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg
  ), list(mean=mean, sd=sd))

sumres[['SumTime']] <- NA
for (m in unique(sumres$Method)) {
  times <- sumres$Time_mean[sumres$Method==m]
  sumres$SumTime[sumres$Method==m] <- cumsum(times)
}

q <- 3
N0 <- 5000 * (1:3)
N1 = 5000


tabres <- ungroup(sumres) |>
  filter(stringr::str_detect(Method, "Batch") | N %in% N0) |>
  dplyr::select(Method, N, SumTime, RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean)
for (i in 1:nrow(tabres)) {
  cat(paste(c(as.character(tabres[i,1:2]), round(as.numeric(tabres[i,3]),1),
    round(as.numeric(tabres[i,4:6]),3)), collapse=" & "), "\\\\ \n")
}

tabres = tabres |> filter(Method != "Batch-PACE")
# PACE is optimized by C codes. Performance too different



# FPC plots --------------------------------------------------------------

source("./R/fdaMdim.R")
source("./data_generation/generator.R")

dat <- gendata.aqi(1)
tgrid <- evalGrid <- dat$tgrid
t0 <- dat$t0
t1 <- dat$t1
neval <- nrow(evalGrid)
muTrueEval <- numeric(neval)
PhiTrueEval <- dat$Phi
lambdaTrue <- dat$lam
q <- npc <- ncol(PhiTrueEval)

basis <- TensorBasis(list(
  create.bspline.basis(c(0, 1), nbasis = 6, norder = 4),
  create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
))
p <- attr(basis, "nbasis")

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
  Method="Truth",
  t1=evalGrid[,1],
  t2=evalGrid[,2],
  phi1=PhiTrueEval[,1],
  phi2=PhiTrueEval[,2],
  phi3=PhiTrueEval[,3]
)
dat.sgd <- data.frame(
  Method="OnlineFPCA-SGD",
  t1=evalGrid[,1],
  t2=evalGrid[,2],
  phi1=PhiSGD[[nPass]][,1],
  phi2=PhiSGD[[nPass]][,2],
  phi3=PhiSGD[[nPass]][,3]
)
dat.adam <- data.frame(
  Method="OnlineFPCA-Adam",
  t1=evalGrid[,1],
  t2=evalGrid[,2],
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
      phi1 = "phi[1](t[1], t[2])",
      phi2 = "phi[2](t[1], t[2])",
      phi3 = "phi[3](t[1], t[2])"
    )
  ) |> 
  mutate(Method =
    recode(Method,
      Truth            = "bold(Truth)",
      `OnlineFPCA-SGD`   = "bold(OnlineFPCA-SGD)",
      `OnlineFPCA-Adam`  = "bold(OnlineFPCA-Adam)"
    )
  )

# Modify the plots to remove legends
ggplot(plotdat) +
  geom_raster(aes(x = t1, y = t2, fill = phiVal), interpolate = TRUE) +
  facet_grid(
    rows     = vars(phiId),
    cols     = vars(Method),
    labeller = label_parsed     # interpret phi[1], phi[2], phi[3] as expressions
  ) +
  theme_bw() +
  theme(
    panel.grid    = element_blank(),
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text       = element_text(colour = "black", face = "bold")
  ) +
  labs(
    x    = expression(t[1]),    # math axis labels
    y    = expression(t[2]),
    fill = NULL                 # drop legend title
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0.0, 0.5, 1.0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0.0, 0.5, 1.0)) +
  scale_fill_viridis_c()        # continuous viridis palette

ggsave(
  file.path(dirpath, "simu-aqi-fpc.pdf"), device = "pdf",
  width = 6, height = 4.8
)
