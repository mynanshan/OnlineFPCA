library(dplyr)
library(tidyr)
library(readr)
library(foreach)

exprmt <- "abv1d"
dirpath <- file.path("experiments", exprmt)


seed <- 1
n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("simu_",exprmt,"_sd",seedtext,".csv")

nBlockList <- c(10, 20, 50, 100)
ewmabv.beta.list <- c(0.1, 0.3, 0.5, 0.7, 0.9)

tau_paths = foreach(nBlock = nBlockList) %:%
  foreach(ewmabv.beta = ewmabv.beta.list, .combine = cbind) %do% {
    spectxt = paste0("nBlock",nBlock,"_beta",ewmabv.beta*10)
    tau_path = readRDS(file.path(dirpath, paste0(
      "taupath_",spectxt,"_",seedtext,".rds"
    )))
    tau_path = matrix(tau_path, ncol=1)
    colnames(tau_path) <- spectxt
    tau_path
  }


N1 <- 5000
N0 <- c(1:3) * N1


res <- read_csv(file.path(dirpath, filename), show_col_types = FALSE) |> 
  filter(N %in% N0) |>
  mutate(
    seed = as.integer(seed),
    epoch = as.integer(N / N1)
  ) |> 
  arrange(N, nBlock, wBeta) |> 
  relocate(epoch, .after = N) |> 
  dplyr::select(epoch, nBlock, wBeta, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg)

do.call(cbind, lapply(tau_paths, tail, n=1))[1,]
