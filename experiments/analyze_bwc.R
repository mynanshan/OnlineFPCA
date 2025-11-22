library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

pa <- par(no.readonly = TRUE)

exprmt <- "bwc"
dirpath <- file.path("experiments", exprmt)

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

sumres <- res %>% group_by(Method, N, C, W) %>%
  summarise_at(vars(
    RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg
  ), list(mean=mean, sd=sd))

# sumres[['SumTime']] <- NA
# for (m in unique(sumres$Method)) {
#   times <- sumres$Time_mean[sumres$Method==m]
#   sumres$SumTime[sumres$Method==m] <- cumsum(times)
# }

q <- 3
N0 <- 5000 * (1:3)
N1 = 5000


# sumres |> 
#   filter(N == N1) |> 
#   group_by(Method, C, W) |> 
#   summarise_at(vars(
#     RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
#   ), list(mean=mean))

# sumres |> 
#   filter(N == 2 * N1) |> 
#   group_by(Method, C, W) |> 
#   summarise_at(vars(
#     RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
#   ), list(mean=mean))

# sumres |> 
#   filter(N == 3 * N1) |> 
#   group_by(Method, C, W) |> 
#   summarise_at(vars(
#     RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
#   ), list(mean=mean))


tabres <- sumres |> 
  filter(N %in% N0, Method == "OnlineFPCA-adagrad") |> 
  group_by(Method, N, C, W) |> 
  summarise_at(vars(
    RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
  ), list(mean=mean)) |> 
  ungroup() |> 
  rename(
    phi1 = RMSEphi1.avg_mean_mean,
    phi2 = RMSEphi2.avg_mean_mean,
    phi3 = RMSEphi3.avg_mean_mean
  ) |> 
  mutate(
    Epoch = N / 5000,
    B = C / W
  ) |> 
  relocate(B, .before = phi1) |> 
  relocate(Epoch, .before = C) |> 
  dplyr::select(-N, -Method) |>
  pivot_wider(
    id_cols = c("C", "W", "B"),
    names_from = "Epoch",
    values_from = c("phi1", "phi2", "phi3")
  )

knitr::kable(tabres, format = "latex", digits = 3)
