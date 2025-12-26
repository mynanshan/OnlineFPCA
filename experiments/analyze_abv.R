library(dplyr)
library(tidyr)
library(readr)

exprmt <- "abv1d"
dirpath <- file.path("experiments", exprmt)


N1 <- 5000
N0 <- c(1:3) * N1


# Check average performance ----------------------------------------------

res <- do.call(
  rbind,
  lapply(
    1:100, \(seed) {
      n_digit_seed = ceiling(log10(seed + 1))
      seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
      filename <- paste0("simu_",exprmt,"_sd",seedtext,".csv")
      read_csv(file.path(dirpath, filename), show_col_types = FALSE)
    }
  )
)

res <- res |> 
  mutate(
    epoch = as.integer(epoch),
    nBlock = as.integer(nBlock)
  ) |> 
  dplyr::select(epoch, nBlock, wBeta, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg, tau) |>
  group_by(epoch, nBlock, wBeta) |> 
  summarise_at(c("RMSEphi1.avg", "RMSEphi2.avg", "RMSEphi3.avg", "tau"), mean) |> 
  arrange(epoch, nBlock, wBeta) 

knitr::kable(res, "latex", digits = 3)

tabres <- res |> ungroup() |> 
  filter(nBlock == 100) |> 
  rename(
    phi1 = RMSEphi1.avg,
    phi2 = RMSEphi2.avg,
    phi3 = RMSEphi3.avg
  ) |> 
  dplyr::select(-nBlock) |>#1
  pivot_wider(
    id_cols = c("wBeta"),
    names_from = "epoch",
    values_from = c("phi1", "phi2", "phi3")
  )

knitr::kable(tabres, format = "latex", digits = 3)
