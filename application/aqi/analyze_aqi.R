library(fda)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/onlineFPCA.R")

### Settings
datapath <- "data/epa-aqs"
respath <- "application"

load(file.path(datapath, "aqi-us.Rda"))
siteinfo <- readxl::read_excel(file.path(datapath, "siteinfo.xlsx"))

lonrange <- c(-125, -67)
latrange <- c(25, 49)
yearrange <- 1982:2022
latgrid <- seq(latrange[1] + 1, latrange[2], 2)
longrid <- seq(lonrange[1] + 1, lonrange[2], 2)
locGrid <- margins2grid(list(latgrid, longrid))
# save(latrange, lonrange, yearrange, locGrid, evalGrid, basis, file="aqi_settings.Rdata")
colnames(locGrid) <- c("lat", "lon")

locGridRescale <- scale(locGrid,
  center = c(latrange[1], lonrange[1]),
  scale = c(diff(latrange), diff(lonrange))
)

evalGridList <- list(
  (latgrid - latrange[1]) / diff(latrange),
  (longrid - lonrange[1]) / diff(lonrange)
)
evalGrid <- margins2grid(evalGridList)

# numbering the locations
lat_id <- match(dat$Latitude.Binned, latgrid)
lon_id <- match(dat$Longitude.Binned, longrid)
dat$LocId <- (lon_id - 1) * length(latgrid) + lat_id

basis <- TensorBasis(list(
  create.bspline.basis(c(0, 1), nbasis = 6, norder = 4),
  create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
))
p <- attr(basis, "nbasis")

N <- length(unique(dat$Date.Local))
q <- 6
nBatch <- 10
nParams <- 6
nPass <- 10
nBlock <- 100
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
# stepsizeList <- c(1e0, 3e-1, 1e-1, 3e-2)
stepsize <- 2e-1
stepsize.min <- 5e-2
sgd.step.scale <- 1 # use a smaller step size for sgd
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
asgd.use <- TRUE
nBlockIter <- nBlock / nBatch
nIter1pass <- round(N / nBatch)
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)


### run the uncommented codes from aqi.R
respath <- "application"
load(file.path(respath, "fit_aqi.Rdata"))
load(file.path(respath, "fit_aqi_batch.Rdata"))

# online results
tau.min <- fit$tau.min
l <- which(fit$tau == tau.min)[1]
Theta.avg <- fit$Theta.avg[, , l]
lambda.avg <- fit$lambda.avg[, l]
sigma2.avg <- fit$sigma2.avg[l]

B <- eval_basis(locGridRescale, basis)
G <- get_basis_inprod_matrix(basis)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)

ord <- order(lambda.avg, decreasing = TRUE)
PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg[, ord], basis))
# approximated FVE
fve <- round(lambda.avg[ord] / sum(lambda.avg), 3) * 100

params <- fit$params.history$params
vcrits <- fit$vcrit.history

# check eigenfunctions
tau_path <- with(fit$tau.select, extract_tau_path(tau.history, tau.selectId, l))
tau_path_id <- c(tau_path$tau_path_id, rep(1, (nPass - 1) * round(N / nBlock)))
tau_path_id_extend <- c(rep(tau_path_id[1], nRoundNoTune), tau_path_id)

pdf(file = file.path(respath, "aqi-taupath.pdf"), width = 7, height = 4.5)
with(fit$tau.select, plot.tau_path2(tau.history, tau.selectId, l))
dev.off()

ThetaAll.avg <- sapply(seq_along(fit$params.history$iter.params),
  \(i) params$Theta.avg[, , tau_path_id_extend[i], i],
  simplify = "array"
)

saveRDS(PhiAvgEval, file.path(respath, "eigf_aqi.rds"))
saveRDS(lambda.avg[ord], file.path(respath, "eigval_aqi.rds"))
saveRDS(evalGrid, file.path(respath, "tgrid_aqi.rds"))

# plot the results ===============
library(ggplot2)
library(sf)

nfpc <- 4

data(us_states, package = "spData")
espg <- as.numeric(stringr::str_extract(st_crs(us_states)$input, "\\d+"))

# create a background geometry
backrect <- list(lonrange + c(-20, 20), latrange + c(-10, 10)) %>%
  margins2grid() %>%
  .[c(1, 2, 4, 3, 1), ] %>%
  list() %>%
  st_polygon() %>%
  st_sfc(crs = espg) %>%
  st_sf()
mapmask <- st_difference(backrect, st_union(us_states))
us_bbox <- st_bbox(st_union(us_states))


plot_fpc_map <- function(Phi, fve) {
  plotdat <- data.frame(
    Latitude = locGrid[, "lat"], Longitude = locGrid[, "lon"]
  )
  for (k in 1:nfpc) {
    eigvals <- Phi[, k]
    plotdat[[paste0("phi", k)]] <- eigvals
  }

  # Modify the plots to remove legends
  figs <- lapply(1:nfpc, \(k) {
    fig <- ggplot() +
      geom_raster(aes_string(x = "Longitude", y = "Latitude", fill = paste0("phi", k)),
        data = plotdat, interpolate = TRUE
      ) +
      geom_sf(data = mapmask, fill = "white") +
      geom_sf(data = us_states, alpha = 0) +
      lims(x = lonrange + c(-1, 1), y = latrange + c(-1, 1)) +
      theme_bw() +
      ggtitle(paste0("FPC ", k, " (", fve[k], "%)")) +
      scale_fill_continuous(
        type = "viridis",
        limits = c(min(Phi), max(Phi))
      ) +
      theme(legend.position = "none") +
      labs(x = NULL, y = NULL)
    return(fig)
  })
  return(figs)
}

figs <- plot_fpc_map(PhiAvgEval[, 1:nfpc], fve)
# figs_batch = plot_fpc_map(PhiBatchEval[,1:nfpc], fve_batch)

# Create a plot for the legend
legend_plot <- ggplot() +
  geom_raster(aes(x = 1, y = 1, fill = 1), show.legend = TRUE) +
  scale_fill_continuous(
    type = "viridis",
    limits = c(min(PhiAvgEval), max(PhiAvgEval))
  ) +
  theme_void() +
  labs(fill = "FPC Value") +
  theme(legend.key.height = unit(2, "cm"))

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

# Combine your plots and the legend
plots_grid <- cowplot::plot_grid(plotlist = figs, ncol = 2, nrow = 2)
# cowplot::plot_grid(plotlist = figs_batch, ncol = 1, nrow = 3)
final_plot <- cowplot::plot_grid(plots_grid, legend, ncol = 2, rel_widths = c(1, 0.2))

final_plot

ggsave(file.path(respath, "aqi-eigfun.pdf"), width = 10.8, height = 5.4)


### running time
# online
times <- c(
  colSums(fit$time.history[, 1:nIter1pass]),
  fit$time.history[1, (nIter1pass + 1):(nPass * nIter1pass)]
)
times <- colSums(matrix(times, nrow = nBlockIter))
times <- c(time_init, times)
sum(times)
# batch
batch_time
