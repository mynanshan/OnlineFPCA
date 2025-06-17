### run the uncommented codes from aqi.R
respath <- "application"
load(file.path(respath, "fit_aqi.Rdata"))
load(file.path(respath, "fit_aqi_batch.Rdata"))

# online results
tau.min <- fit$tau.min
l <- which(fit$tau == tau.min)[1]
Theta.avg <- fit$Theta.avg[,,l]
lambda.avg <- fit$lambda.avg[,l]
sigma2.avg <- fit$sigma2.avg[l]

ord = order(lambda.avg, decreasing = TRUE)
PhiAvgEval <- eval_fd(evalGrid, FuncData(Theta.avg[,ord], basis))
# approximated FVE
fve = round(lambda.avg[ord] / sum(lambda.avg), 3) * 100

params <- fit$params.history$params
vcrits <- fit$vcrit.history

# check eigenfunctions
tau_path <- with(fit$tau.select, extract_tau_path(tau.history, tau.selectId, l))
tau_path_id <- c(tau_path$tau_path_id, rep(1, (nPass-1)*round(N/nBlock)))
tau_path_id_extend <- c(rep(tau_path_id[1], nRoundNoTune), tau_path_id)

setEPS()
postscript(file = file.path(respath, "aqi-taupath.eps"), width=7, height=4.5)
with(fit$tau.select, plot.tau_path2(tau.history, tau.selectId, l))
dev.off()

ThetaAll.avg <- sapply(seq_along(fit$params.history$iter.params),
  \(i) params$Theta.avg[,,tau_path_id_extend[i],i],
  simplify = "array")
# PhiAvgAll <- eval_fd(evalGrid, FuncData(ThetaAll.avg, basis))

# # batch results
# ThetaBatch <- fitBatch$Theta
# muBatchEval <- eval_fd(evalGrid, FuncData(fitBatch$theta_mu, basis))
# PhiBatchEval <- eval_fd(evalGrid, FuncData(ThetaBatch, basis))
# lambdaBatch = colMeans(fitBatch$alphamat^2)
# fve_batch = round(lambdaBatch / sum(lambdaBatch), 3) * 100


# plot the results ===============
library(ggplot2)
library(sf)

nfpc = 4

data(us_states, package = "spData")
espg <- as.numeric(stringr::str_extract(st_crs(us_states)$input, "\\d+"))

# create a background geometry
backrect = list(lonrange+c(-20,20), latrange+c(-10,10)) %>%
  margins2grid() %>% .[c(1,2,4,3,1),] %>%
  list() %>% st_polygon() %>%
  st_sfc(crs=espg) %>% st_sf()
mapmask = st_difference(backrect, st_union(us_states))
us_bbox = st_bbox(st_union(us_states))


plot_fpc_map = function(Phi, fve) {
  
  plotdat <- data.frame(
    Latitude=locGrid[,"lat"], Longitude=locGrid[,"lon"]
  )
  for (k in 1:nfpc) {
    eigvals = Phi[,k]
    plotdat[[paste0('phi',k)]] = eigvals
  }
  
  # Modify the plots to remove legends
  figs = lapply(1:nfpc, \(k) {
    fig = ggplot() +
      geom_raster(aes_string(x="Longitude", y="Latitude", fill=paste0("phi",k)),
                  data=plotdat, interpolate = TRUE) +
      geom_sf(data=mapmask, fill="white") +
      geom_sf(data=us_states, alpha=0) +
      lims(x=lonrange+c(-1,1), y=latrange+c(-1,1)) +
      theme_bw() +
      ggtitle(paste0("FPC ", k, ", ", fve[k], "%")) +
      scale_fill_continuous(type = "viridis",
                            limits=c(min(Phi), max(Phi))) +
      theme(legend.position = "none") +
      labs(x = NULL, y = NULL)
    return(fig)
  })
  return(figs)
}

figs = plot_fpc_map(PhiAvgEval[,1:nfpc], fve)
# figs_batch = plot_fpc_map(PhiBatchEval[,1:nfpc], fve_batch)

# Create a plot for the legend
legend_plot <- ggplot() +
  geom_raster(aes(x=1, y=1, fill=1), show.legend = TRUE) +
  scale_fill_continuous(type = "viridis",
    limits=c(min(PhiAvgEval), max(PhiAvgEval))) +
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

ggsave(file.path(respath, "aqi-eigfun.eps"), width=10.8, height=5.4)


### running time
# online
times <- c(colSums(fit$time.history[,1:nIter1pass]),
           fit$time.history[1,(nIter1pass+1):(nPass*nIter1pass)])
times <- colSums(matrix(times, nrow=nBlockIter))
times <- c(time_init, times)
sum(times)
# batch
batch_time

# AQI surface reconstruction ==============================

# tmpdat <- dat |> filter(subjId==1503)
# X <- tmpdat |> dplyr::select(s, t) |> as.matrix()
# y <- tmpdat |> pull(z)
# 
# 
# Phi <- eval_fd(X, FuncData(ThetaInit, basis))
# # fit.subj <- glmnet::cv.glmnet(Phi, y, penalty.factor = 1/lambdaInit, intercept=FALSE)
# # lambda.min <- fit.subj$lambda.min
# fit.subj <- glmnet::glmnet(Phi, y, penalty.factor = 1/lambdaInit, lambda = 1e-6, intercept=FALSE)
# fit.subj$beta
# tmpdat$zpred <- predict(fit.subj, newx = Phi)
# # fit.subj <- lm(y~Phi-1)
# # tmpdat$zpred <- predict(fit.subj)
# 
# ggplot() +
#   geom_point(aes(x=Longitude.Binned, y=Latitude.Binned,
#     color=z, size=z), stroke=1,
#     data=tmpdat, shape=1) +
#   geom_sf(data=mapmask, fill="white")
# geom_sf(data=us_states) +
#   theme_bw() +
#   scale_color_continuous(type = "viridis") +
#   labs(color = "LogAQI", size = "LogAQI")
# 
# ggplot() +
#   geom_sf(data=us_states) +
#   geom_point(aes(x=Longitude.Binned, y=Latitude.Binned,
#     color=zpred, size=zpred), stroke=1,
#     data=tmpdat, shape=1) +
#   theme_bw() +
#   scale_color_continuous(type = "viridis") +
#   labs(color = "LogAQI", size = "LogAQI")
# 
# 
# load(file.path(respath, "fit_aqi.Rdata"))
# colnames(locGrid) = c("lat", "lon")
# 
# plotdat <- data.frame(
#   Latitude=locGrid[,"lat"], Longitude=locGrid[,"lon"],
#   InUS=in_us_idx
# )
# for (k in 1:4) {
#   eigvals = PhiAvgEval[,k]
#   # eigvals[!in_us_idx] = NA
#   plotdat[[paste0('phi',k)]] = eigvals
# }
# 
# figs = lapply(1:4, \(k) {
#   fig = ggplot() +
#     geom_raster(aes_string(x="Longitude", y="Latitude", fill=paste0("phi",k)),
#       data=plotdat, interpolate = TRUE) +
#     geom_sf(data=mapmask, fill="white") +
#     geom_sf(data=us_states, alpha=0) +
#     lims(x=lonrange+c(-1,1), y=latrange+c(-1,1)) +
#     theme_bw() +
#     ggtitle(paste("FPC", k)) +
#     scale_fill_continuous(type = "viridis",
#       limits=c(min(PhiAvgEval), max(PhiAvgEval))) +
#     theme(legend.title=element_blank())
#   return(fig)
# })
# 
# cowplot::plot_grid(plotlist = figs,nrow=2,ncol=2)

