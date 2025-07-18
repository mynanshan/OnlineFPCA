library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

exprmt <- "perm1d"
dirpath <- file.path("experiments", exprmt)


# RMSE uncertainty -------------------------------------------------------


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

q <- 3
N0 <- 5000 * (1:3)
N1 = 5000

plotdat <- res |> 
  mutate(n = N) |> 
  # mutate(n = as.factor(N)) |> 
  pivot_longer(
    cols = paste0("RMSEphi", 1:q, ".avg"),
    names_to = "phiId",
    values_to = "RMSE"
  ) |> 
  mutate(
    phiId = recode(
      phiId,
      RMSEphi1.avg = "phi[1](t)",
      RMSEphi2.avg = "phi[2](t)",
      RMSEphi3.avg = "phi[3](t)"
    )
  )


ggplot(plotdat) +
  geom_boxplot(
    aes(x = n, y = RMSE, group=n),
    linewidth     = 0.3,  # thinner box outlines (default is 0.5)
    outlier.size  = 0.2   # smaller outlier points (default is 2.5)
  ) +
  facet_wrap(
    vars(phiId),
    labeller = label_parsed,
    nrow = 3,
    scales = "free_y",
    strip.position = "right"
  ) +
  # scale_y_log10() +
  theme_bw() +
  theme(
    strip.placement = "outside"     # put the strip outside the panel area
  )


# ggsave(file.path(dirpath, "perm1d-rmse.pdf"), device="pdf",
#   width = 7, height = 6)


# FPC UQ Band Plot -------------------------------------------------------

source("./R/fdaMdim.R")
source("./data_generation/generator.R")

alpha = 2
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
basis <- create.bspline.basis(c(t0,t1), nbasis = nbasis, norder = 4)
p <- nbasis

file_list = list.files(dirpath)
fpc_files = file_list[startsWith(file_list, "Theta")]
nrep = length(fpc_files)

nPass <- 3
ThetaSGD <- ThetaAdam <- rep(list(array(0, c(p, q, nrep))), 3)
for (irep in seq_along(fpc_files)) {
  ThetaRecord <- readRDS(file.path(dirpath, fpc_files[irep]))  
  for (ip in seq_len(nPass)) {
    ThetaSGD[[ip]][,,irep] <- ThetaRecord[[ip]][['SGD']]
    ThetaAdam[[ip]][,,irep] <- ThetaRecord[[ip]][['Adam']]
  }
}

PhiSGD <- lapply(ThetaSGD, \(Theta) {
  apply(
    Theta, 3,
    \(X) eval_fd(evalGrid, FuncData(X, basis))
  ) |> # (nt*q, nrep)
    array(c(neval, q, nrep))
})
PhiAdam <- lapply(ThetaAdam, \(Theta) {
  apply(
    Theta, 3,
    \(X) eval_fd(evalGrid, FuncData(X, basis))
  ) |> # (nt*q, nrep)
    array(c(neval, q, nrep))
})


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
    seq_len(nrep), \(irep) {
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
dat.adam <- do.call(
  bind_rows,
  lapply(
    seq_len(nrep), \(irep) {
      data.frame(
        irep=irep,
        Method="OnlineFPCA-Adam",
        t=evalGrid,
        phi1=PhiAdam[[nPass]][,1,irep],
        phi2=PhiAdam[[nPass]][,2,irep],
        phi3=PhiAdam[[nPass]][,3,irep]
      )
    }
  )
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


# 1) pick off the truth curves
truth_curve <- plotdat %>%
  filter(Method == "Truth") %>%
  select(-Method)

# 2) replicate for each method‐row you want
methods_to_show <- c("OnlineFPCA-SGD", "OnlineFPCA-Adam")
truth_for_all <- crossing(truth_curve, Method = methods_to_show) |> 
  mutate(Method = factor(Method, levels = c("OnlineFPCA-SGD", "OnlineFPCA-Adam")))

# 3) build the plot
ggplot() +
  # the cloud of replicates, only for the two methods
  geom_line(
    data  = plotdat %>% filter(Method != "Truth") |> 
      mutate(Method = factor(Method, levels = c("OnlineFPCA-SGD", "OnlineFPCA-Adam"))),
    aes(x = t, y = phiVal, group = irep),
    alpha = 0.1, color = "#81a9d3",
    linewidth = 2
  ) +
  # the red truth curve, now tagged onto each method level
  geom_line(
    data  = truth_for_all,
    aes(x = t, y = phiVal, group = irep),
    color = "#cb3064",
    # color = "black",
    # alpha = 0.8,
    linewidth = 0.8
  ) +
  facet_grid(
    rows     = vars(Method),
    cols     = vars(phiId),
    labeller = label_parsed
  ) +
  theme_bw() +
  theme(
    panel.grid    = element_blank(),
    axis.text.x   = element_text(size = 8)
  ) +
  labs(
    x = expression(t),
    y = "FPC Value"
  )


ggsave(
  file.path(dirpath, "perm1d-fpc.pdf"), device = "pdf",
  width = 6, height = 4
)
