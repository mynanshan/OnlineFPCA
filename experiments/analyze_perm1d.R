library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

exprmt <- "perm1d"
dirpath <- file.path("experiments", exprmt)


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

nrec <- 3 * 5
ndata.1rec <- 5000 / 5
ThetaSGD <- ThetaAdam <- rep(list(array(0, c(p, q, nrep))), nrec)
for (irep in seq_along(fpc_files)) {
  ThetaRecord <- readRDS(file.path(dirpath, fpc_files[irep]))  
  for (ip in seq_len(nrec)) {
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
  n=0,
  Method="Truth",
  t=evalGrid,
  phi1=PhiTrueEval[,1],
  phi2=PhiTrueEval[,2],
  phi3=PhiTrueEval[,3]
)
dat.adam <- do.call(
  bind_rows,
  lapply(
    seq_len(nrec), \(irec) {
      do.call(
        bind_rows,
        lapply(
          seq_len(nrep), \(irep) {
            data.frame(
              irep=irep,
              n=irec * ndata.1rec,
              Method="OnlineFPCA-Adam",
              t=evalGrid,
              phi1=PhiAdam[[irec]][,1,irep],
              phi2=PhiAdam[[irec]][,2,irep],
              phi3=PhiAdam[[irec]][,3,irep]
            )
          }
        )
      )
    }
  )
)
plotdat <- bind_rows(dat.true, dat.adam) %>%
  pivot_longer(
    cols = starts_with("phi"),
    names_to  = "phiId",
    values_to = "phiVal"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("Truth", "OnlineFPCA-Adam")),
    # map “phi1”→“phi[1]” etc so label_parsed knows to render φ₁, φ₂, φ₃
    phiId = recode(phiId,
      phi1 = "phi[1](t)",
      phi2 = "phi[2](t)",
      phi3 = "phi[3](t)"
    )
  )

# ndata.levels <- unique(plotdat$n)
ndata.levels <- c(1000, 5000, 10000)


# 1) pick off the truth curves
truth_curve <- plotdat %>%
  filter(Method == "Truth") %>%
  dplyr::select(-Method, -n)

# 2) replicate for each method‐row you want
truth_for_all <- crossing(truth_curve, n = ndata.levels)

my_n_labeller <- function(x) {
  # x comes in as character/factor; force it numeric
  num_x <- as.numeric(as.character(x))
  # format with commas and prefix “n=”
  paste0("n=", comma(num_x))
}

# 3) build the plot
ggplot() +
  # the cloud of replicates, only for the two methods
  geom_line(
    data  = plotdat %>% filter(Method != "Truth", n %in% ndata.levels),
    aes(x = t, y = phiVal, group = irep),
    alpha = 0.1, color = "#81a9d3",
    linewidth = 1
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
    rows     = vars(n),
    cols     = vars(phiId),
    labeller = labeller(
      # for the n‑strips: format with commas and prefix “n=”
      n      = my_n_labeller,
      # for the phi‑strips: parse expressions like "phi[1]" → φ₁
      phiId  = label_parsed
    )
  ) +
  theme_bw() +
  theme(
    panel.grid    = element_blank(),
    axis.text.x   = element_text(size = 8),
    strip.background = element_blank()
  ) +
  labs(
    x = expression(t),
    y = "FPC Value"
  )


ggsave(
  file.path(dirpath, "perm1d-fpc.pdf"), device = "pdf",
  width = 6, height = 4
)
