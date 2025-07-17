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
  geom_boxplot(aes(x = n, y = RMSE, group=n)) +
  facet_wrap(
    vars(phiId),
    labeller = label_parsed,
    nrow = 3,
    scales = "free_y"
  ) +
  scale_y_log10() +
  theme_bw()
