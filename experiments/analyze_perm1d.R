
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

pa <- par(no.readonly = TRUE)

exprmt <- "perm1d"
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

q <- 3
N0 <- 5000 * (1:3)
N1 = 5000


ggplot(res |> mutate(N == as.factor(N))) +
  geom_boxplot(aes(x = N, y = RMSEphi1.avg, group=N))

ggplot(res |> mutate(N == as.factor(N))) +
  geom_boxplot(aes(x = N, y = RMSEphi2.avg, group=N))

ggplot(res |> mutate(N == as.factor(N))) +
  geom_boxplot(aes(x = N, y = RMSEphi3.avg, group=N))
