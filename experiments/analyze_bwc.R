library(tidyr)
library(dplyr)
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

sumres <- res %>% group_by(Method, N, B, W) %>%
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


sumres |> 
  filter(N == N1) |> 
  group_by(Method, B, W) |> 
  summarise_at(vars(
    RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
  ), list(mean=mean, sd=sd))

sumres |> 
  filter(N == 2 * N1) |> 
  group_by(Method, B, W) |> 
  summarise_at(vars(
    RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
  ), list(mean=mean, sd=sd))

sumres |> 
  filter(N == 3 * N1) |> 
  group_by(Method, B, W) |> 
  summarise_at(vars(
    RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
  ), list(mean=mean, sd=sd))
