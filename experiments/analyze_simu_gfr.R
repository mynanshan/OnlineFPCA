library(dplyr)
library(tidyr)
library(readr)
library(fda)
library(ggplot2)

exprmt <- "gfr"
dirpath <- file.path("experiments", exprmt)


# RMSE plot and table ----------------------------------------------------

res <- do.call(
  rbind,
  lapply(
    1:100, \(i) {
      seed <- i
      n_digit_seed <- ceiling(log10(seed + 1))
      seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
      filename <- paste0("simu_", exprmt, "_sd", seedtext, ".csv")
      tryCatch(
        {
          read_csv(file.path(dirpath, filename), show_col_types = FALSE)
        },
        error = function(e) {
          NULL
        }
      )
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

sumres <- res |>
  group_by(Method, StepSize, Ninit, npc, initMethod, N) |>
  summarise_at(vars(Time, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg), list(mean = mean)) |>
  ungroup() |>
  group_by(Method, StepSize, Ninit, npc, initMethod) |>
  mutate(SumTime = cumsum(Time_mean))

q <- 3
N0 <- c(5000, 10000, 15000)
N1 <- 5000

meth_ord <- c(
  paste0("OnlineFPCA-", c("sgd", "adagrad")),
  "OnlineCov",
  "Batch-PACE", "Batch-FACE", "Batch-REML", "Batch-SOAP"
)

tabres <- ungroup(sumres) |>
  filter(stringr::str_detect(Method, "Batch") | N %in% N0) |>
  dplyr::select(
    Method, N, StepSize, SumTime,
    RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
  ) |>
  mutate(Method = factor(Method, levels = meth_ord)) |>
  arrange(Method, N) |>
  mutate(
    epoch = as.integer(N / N1),
    Method = case_match(
      Method,
      "OnlineFPCA-sgd" ~ "OnlineFPCA-RSGD",
      "OnlineFPCA-adagrad" ~ "OnlineFPCA-RAdagrad",
      "OnlineCov" ~ "OnlineCov",
      "Batch-PACE" ~ "PACE",
      "Batch-FACE" ~ "FACE",
      "Batch-REML" ~ "REML",
      "Batch-SOAP" ~ "SOAP"
    ),
    SumTime = round(SumTime, 1)
  ) |>
  rename(
    RMSEphi1 = RMSEphi1.avg_mean,
    RMSEphi2 = RMSEphi2.avg_mean,
    RMSEphi3 = RMSEphi3.avg_mean
  ) |>
  relocate(epoch, .before = SumTime) |>
  dplyr::select(-N)

knitr::kable(
  tabres |>
    filter(is.na(StepSize) | StepSize == 0.1) |>
    dplyr::select(-StepSize),
  "latex",
  digits = 3
)


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
basis <- create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = 4)
p <- nbasis

file_list <- list.files(dirpath)
fpc_files <- file_list[startsWith(file_list, "Theta")]
nrep <- length(fpc_files)

nPass <- 3
ThetaSGD <- ThetaAdaGrad <- rep(list(matrix(0, nrow = p, ncol = q)), 3)
for (f in fpc_files) {
  tmp <- tryCatch(
    readRDS(file.path(dirpath, f)),
    error = function(e) {
      NULL
    }
  )
  if (is.null(tmp)) next
  ThetaRecord <- tmp
  for (ip in seq_len(nPass)) {
    ThetaSGD[[ip]] <- ThetaSGD[[ip]] + ThetaRecord[[ip]][["sgd"]] / nrep
    ThetaAdaGrad[[ip]] <- ThetaAdaGrad[[ip]] + ThetaRecord[[ip]][["adagrad"]] / nrep
  }
}

PhiSGD <- lapply(ThetaSGD, \(Theta) eval_fd(evalGrid, FuncData(Theta, basis)))
PhiAdaGrad <- lapply(ThetaAdaGrad, \(Theta) eval_fd(evalGrid, FuncData(Theta, basis)))


dat.true <- data.frame(
  irep = 0,
  Method = "Truth",
  t = evalGrid,
  phi1 = PhiTrueEval[, 1],
  phi2 = PhiTrueEval[, 2],
  phi3 = PhiTrueEval[, 3]
)
dat.sgd <- data.frame(
  Method = "OnlineFPCA-RSGD",
  t = evalGrid,
  phi1 = PhiSGD[[nPass]][, 1],
  phi2 = PhiSGD[[nPass]][, 2],
  phi3 = PhiSGD[[nPass]][, 3]
)
dat.adagrad <- data.frame(
  Method = "OnlineFPCA-RAdaGrad",
  t = evalGrid,
  phi1 = PhiAdaGrad[[nPass]][, 1],
  phi2 = PhiAdaGrad[[nPass]][, 2],
  phi3 = PhiAdaGrad[[nPass]][, 3]
)
plotdat <- bind_rows(dat.true, dat.sgd, dat.adagrad) %>%
  pivot_longer(
    cols = starts_with("phi"),
    names_to = "phiId",
    values_to = "phiVal"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("Truth", "OnlineFPCA-RSGD", "OnlineFPCA-RAdaGrad")),
    # map “phi1”→“phi[1]” etc so label_parsed knows to render φ₁, φ₂, φ₃
    phiId = recode(phiId,
      phi1 = "phi[1](t)",
      phi2 = "phi[2](t)",
      phi3 = "phi[3](t)"
    )
  ) |>
  mutate(
    Method =
      recode(Method,
        Truth = "bold(Truth)",
        `OnlineFPCA-RSGD` = "bold(RSGD)",
        `OnlineFPCA-RAdaGrad` = "bold(RAdaGrad)"
      )
  )


ggplot(plotdat) +
  geom_line(aes(x = t, y = phiVal)) +
  facet_grid(
    rows = vars(phiId),
    cols = vars(Method),
    labeller = label_parsed, # interpret phi[1], phi[2], phi[3] as expressions
    scales = "free_y"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8),
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text = element_text(colour = "black", face = "bold")
  ) +
  labs(
    x    = expression(t), # math axis labels
    y    = expression(bold(FPC)),
    fill = NULL # drop legend title
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0.0, 0.5, 1.0)) +
  scale_fill_viridis_c() # continuous viridis palette

ggsave(
  file.path(dirpath, "simu-gfr-fpc.pdf"),
  device = "pdf",
  width = 5.4, height = 4.8
)
