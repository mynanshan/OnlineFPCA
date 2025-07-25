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
  ) |> 
  mutate(Method =
    recode(Method,
      Truth            = "bold(Truth)",
      `OnlineFPCA-SGD`   = "bold(OnlineFPCA-SGD)",
      `OnlineFPCA-Adam`  = "bold(OnlineFPCA-Adam)"
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
    axis.text.x = element_text(size=8),
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text       = element_text(colour = "black", face = "bold")
  ) +
  labs(
    x    = expression(t),    # math axis labels
    y    = expression(bold(FPC)),
    fill = NULL                 # drop legend title
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0.0, 0.5, 1.0)) +
  scale_fill_viridis_c()        # continuous viridis palette

ggsave(
  file.path(dirpath, "simu-gfr-fpc.pdf"), device = "pdf",
  width = 5.4, height = 4.8
)
