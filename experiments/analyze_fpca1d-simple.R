library(dplyr)
library(tidyr)
library(readr)

exprmt <- "fpca1d-simple"
dirpath <- file.path("experiments", exprmt)


# RMSE table ----------------------------------------------------

noise_sd <- 0.5
sgdtype <- "sgd"

res <- do.call(
  rbind,
  lapply(
    1:100, \(i) {
      seed <- i
      n_digit_seed = ceiling(log10(seed + 1))
      seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
      filename <- paste0(
        "simu_", exprmt, "_sig", noise_sd, "_", sgdtype,
        "_sd", seedtext, ".csv"
      )
      tryCatch(
        {read_csv(file.path(dirpath, filename), show_col_types = FALSE)},
        error = function(e) NULL
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
  # filter(StepSize == 0.1) |> 
  group_by(Method, noise, Ninit, npc, initMethod, N) |> 
  summarise_at(vars(Time, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg), list(mean=mean)) |> 
  ungroup() |> 
  group_by(Method, noise, Ninit, npc, initMethod) |> 
  mutate(SumTime = cumsum(Time_mean))

q <- 3
N0 <- c(5000, 10000, 15000)
N1 = 5000

res |> filter(
  N %in% N0,
  noise == 0.5,
  Method == "OnlineFPCA-sgd"
) |> dplyr::select(
  noise, seed, Method, StepSize, N,
  RMSEphi1, RMSEphi2, RMSEphi3,
  RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg
) |> View()

# PACE is optimized by C codes. Performance too different
meth_ord = c("OnlineFPCA-sgd", "OnlineFPCA-adagrad", "OnlineFPCA-adam", "OnlineFPCA-rasa")

tabres <- ungroup(sumres) |>
  filter(stringr::str_detect(Method, "Batch") | N %in% N0) |>
  dplyr::select(
    Method, N, noise, SumTime,
    RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
  ) |>
  filter(Method != "Batch-PACE") |> 
  mutate(Method = factor(Method, levels=meth_ord)) |> 
  arrange(noise, Method, N) |> 
  relocate(noise, .before = Method) |> 
  mutate(
    epoch = as.integer(N / N1),
    SumTime = round(SumTime, 1)
  ) |> 
  relocate(epoch, .before = SumTime) |> 
  dplyr::select(-N)

knitr::kable(tabres, "latex", digits = 3)
