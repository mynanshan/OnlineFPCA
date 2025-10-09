library(dplyr)
library(tidyr)
library(readr)

exprmt <- "fpca2d"
dirpath <- file.path("experiments", exprmt)


# RMSE table ----------------------------------------------------

res <- do.call(
  rbind,
  lapply(
    1:100, \(i) {
      seed <- i
      n_digit_seed = ceiling(log10(seed + 1))
      seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
      filename <- paste0("simu_",exprmt,"_sd",seedtext,".csv")
      tryCatch(
        read_csv(file.path(dirpath, filename), show_col_types = FALSE),
        error = function(e) NULL
      )
      
    }
  )
)

res <- res |> 
  mutate(
    seed = as.integer(seed),
    nBatch = as.integer(nBatch)
  )

sumres <- res |> 
  group_by(Method, noise, N) |> 
  summarise_at(vars(Time, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg), list(mean=mean)) |> 
  ungroup() |> 
  group_by(Method, noise) |> 
  mutate(SumTime = cumsum(Time_mean))

q <- 5
N1 = 5000
N0 <- N1 * (1:q)

# PACE is optimized by C codes. Performance too different
meth_ord = c("OnlineFPCA-sgd", "OnlineFPCA-adagrad", "Batch-SOAP")

tabres <- ungroup(sumres) |>
  filter(stringr::str_detect(Method, "Batch") | N %in% N0) |>
  dplyr::select(
    Method, N, noise, SumTime,
    RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
  ) |>
  mutate(Method = factor(Method, levels=meth_ord)) |> 
  arrange(noise, Method, N) |> 
  relocate(noise, .before = Method) |> 
  mutate(
    epoch = as.integer(N / N1),
    Method = case_match(Method,
      "OnlineFPCA-sgd" ~ "OnlineFPCA-RSGD",
      "OnlineFPCA-adagrad" ~ "OnlineFPCA-RAdagrad",
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

knitr::kable(tabres, "latex", digits = 3)
