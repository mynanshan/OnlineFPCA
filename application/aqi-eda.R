library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(lubridate)
library(giscoR)
library(sf)
library(spData)

library(fda)
library(Matrix)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./R/onlineFPCA.R")

par.old <- par(no.readonly = TRUE)


data(us_states, package = "spData")
datapath <- "data/epa-aqs"

# lonrange <- c(-125, -60)
# latrange <- c(20, 60)
lonrange <- c(-125, -67)
latrange <- c(25, 49)


yearrange <- 1982:2022
filenameprefix <- "daily_81102_"

usa <- giscoR::gisco_get_countries(country = "USA") |>
  sf::st_crop(c(
    xmin = lonrange[1], xmax = lonrange[2],
    ymin = latrange[1], ymax = latrange[2]
  ))
espg <- as.numeric(stringr::str_extract(st_crs(us_states)$input, "\\d+"))

# variable names
# [1] "State.Code"          "County.Code"         "Site.Num"
# [4] "Parameter.Code"      "POC"                 "Latitude"
# [7] "Longitude"           "Datum"               "Parameter.Name"
# [10] "Sample.Duration"     "Pollutant.Standard"  "Date.Local"
# [13] "Units.of.Measure"    "Event.Type"          "Observation.Count"
# [16] "Observation.Percent" "Arithmetic.Mean"     "X1st.Max.Value"
# [19] "X1st.Max.Hour"       "AQI"                 "Method.Code"
# [22] "Method.Name"         "Local.Site.Name"     "Address"
# [25] "State.Name"          "County.Name"         "City.Name"
# [28] "CBSA.Name"           "Date.of.Last.Change"


# helper functions
crop_latlon <- function(dat) {
  dat |> filter(Longitude > lonrange[1] &
    Longitude < lonrange[2] &
    Latitude > latrange[1] &
    Latitude < latrange[2])
}

bin_latlon <- function(dat, width = 2, suffix = ".Binned",
                       lat.min = -90, lon.min = -180) {
  if (360 %% width != 0) stop("Inappropriate grid width.")
  dat[[paste0("Latitude", suffix)]] <-
    (dat$Latitude - lat.min) %/% width * width + lat.min + width / 2
  dat[[paste0("Longitude", suffix)]] <-
    (dat$Longitude - lon.min) %/% width * width + lon.min + width / 2
  return(dat)
}


## check different sites
siteinfo <- do.call(
  rbind, lapply(file.path(datapath, paste0(filenameprefix, yearrange, ".csv")), \(f) {
    print(f)
    dat <- read.csv(f) |>
      group_by(State.Code, County.Code, Site.Num) |>
      summarise(
        Latitude = unique(Latitude),
        Longitude = unique(Longitude),
        Nrecord = n(), .groups = "drop"
      )
  })
) |>
  group_by(State.Code, County.Code, Site.Num) |>
  summarise(
    Latitude = unique(Latitude),
    Longitude = unique(Longitude),
    Nrecord = sum(Nrecord), .groups = "drop"
  )

siteinfo <- crop_latlon(siteinfo)

ggplot() +
  geom_sf(data = us_states) +
  geom_sf(data = st_sfc(st_point(c(-100, 40)), crs = espg))

st_intersects(us_states, st_sfc(st_point(c(-100, 40)), crs = espg))

ggplot() +
  geom_sf(data = us_states) +
  geom_point(
    aes(
      x = Longitude, y = Latitude, color = Nrecord,
      size = Nrecord, group = Nrecord
    ),
    stroke = 1,
    data = crop_latlon(siteinfo) |> arrange(Nrecord), shape = 1
  ) +
  theme_bw() +
  scale_color_continuous(type = "viridis") +
  labs(color = "Num Records", size = "Num Records")


## Binning the latitudes and longitudes ==========
gridwidth <- 2
rawgrid <- expand.grid(
  list(
    Latitude = seq(latrange[1] + gridwidth / 2, latrange[2] - gridwidth / 2, gridwidth),
    Longitude = seq(lonrange[1] + gridwidth / 2, lonrange[2] - gridwidth / 2, gridwidth)
  )
)
rawgrid[["InUS"]] <- sapply(1:nrow(rawgrid), \(i) {
  ll <- c(rawgrid$Longitude[i], rawgrid$Latitude[i]) |>
    st_point() |>
    st_sfc(crs = espg) |>
    st_intersects(us_states) |>
    unlist() |>
    length()
  return(ll > 0)
})
latlongrid <- rawgrid |>
  filter(InUS) |>
  dplyr::select(Latitude, Longitude)

siteinfo <- siteinfo |>
  bin_latlon(width = gridwidth, lat.min = latrange[1], lon.min = lonrange[1])
siteinfo <- siteinfo |>
  left_join(
    rawgrid |>
      rename(Latitude.Binned = Latitude, Longitude.Binned = Longitude),
    by = c("Latitude.Binned", "Longitude.Binned")
  )
siteinfo[["Latitude.Projected"]] <- siteinfo$Latitude.Binned
siteinfo[["Longitude.Projected"]] <- siteinfo$Longitude.Binned

for (i in 1:nrow(siteinfo)) {
  if (siteinfo$InUS[i]) next
  dist_vec <- (latlongrid$Latitude - siteinfo$Latitude.Binned[i])^2 +
    (latlongrid$Longitude - siteinfo$Longitude.Binned[i])^2
  idx <- which.min(dist_vec)[1]
  siteinfo$Latitude.Projected[i] <- latlongrid$Latitude[idx]
  siteinfo$Longitude.Projected[i] <- latlongrid$Longitude[idx]
}

siteinfo <- siteinfo |>
  mutate(Latitude.Binned = NULL, Longitude.Binned = NULL) |>
  rename(Latitude.Binned = Latitude.Projected, Longitude.Binned = Longitude.Projected)

nrow(siteinfo) # 3443
nrow(siteinfo |> distinct(Latitude.Binned, Longitude.Binned)) # 193

ggplot() +
  geom_sf(data = us_states) +
  geom_point(
    aes(
      x = Longitude.Binned, y = Latitude.Binned, color = Nrecord,
      size = Nrecord, group = Nrecord
    ),
    stroke = 2,
    data = crop_latlon(siteinfo) |>
      group_by(Latitude.Binned, Longitude.Binned) |>
      summarise(Nrecord = sum(Nrecord) / 1000, .groups = "drop") |>
      arrange(Nrecord) |>
      filter(Nrecord > 0),
    shape = 1
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  # scale_color_continuous(type = "viridis") +
  labs(color = "Number of Records (Thousands)", x = NULL, y = NULL) +
  guides(size = "none") +
  theme(legend.key.width = unit(1, "in"))

openxlsx::write.xlsx(siteinfo, file.path(datapath, "siteinfo.xlsx"))


siteinfo <- readxl::read_excel(file.path(datapath, "siteinfo.xlsx"))


## Check num obs per subject

## Check acf
dat <- read.csv(file.path(datapath, "daily_81102_2022.csv"))
dat <- dat |>
  crop_latlon() |>
  bin_latlon() |>
  select(Latitude.Binned, Longitude.Binned, Date.Local, AQI) |>
  group_by(Latitude.Binned, Longitude.Binned, Date.Local) |>
  summarise(AQI = mean(AQI), .groups = "drop") |>
  mutate(Date.Local = ymd(Date.Local)) |>
  arrange(Date.Local)
# TODO: will it be more reasonable to do geometric mean of AQI,
# i.e., arithmetic mean of log AQI?

dat |>
  group_by(Date.Local) |>
  summarise(Nobs = n())

ggplot() +
  geom_sf(data = us_states) +
  geom_point(aes(x = Longitude.Binned, y = Latitude.Binned, color = log1p(AQI)),
    stroke = 1,
    data = dat |> filter(Date.Local == ymd("2022-03-01")), shape = 1
  ) +
  theme_bw() +
  scale_color_continuous(type = "viridis") +
  labs(color = "AQI")


x <- dat |>
  filter(Latitude.Binned == 33, Longitude.Binned == -87) |>
  pull(AQI) |>
  log1p()
acf(diff(x))

y <- x[seq(1, length(x), 2)]
acf(y)

acfcheck <- dat |>
  mutate(logAQI = log1p(AQI)) |>
  group_by(Latitude.Binned, Longitude.Binned) |>
  summarise(
    acf1.max = max(abs(acf(logAQI, plot = F)$acf[-1])),
    acf2.max = max(abs(acf(logAQI[seq(1, n(), 2)], plot = F)$acf[-1])),
    acf5.max = max(abs(acf(logAQI[seq(1, n(), 5)], plot = F)$acf[-1])),
    acf7.max = max(abs(acf(logAQI[seq(1, n(), 7)], plot = F)$acf[-1])),
    acf10.max = max(abs(acf(logAQI[seq(1, n(), 10)], plot = F)$acf[-1])),
    acf1.mean3 = mean(abs(acf(logAQI, plot = F)$acf[2:4])),
    acf2.mean3 = mean(abs(acf(logAQI[seq(1, n(), 2)], plot = F)$acf[2:4])),
    acf5.mean3 = mean(abs(acf(logAQI[seq(1, n(), 5)], plot = F)$acf[2:4])),
    acf7.mean3 = mean(abs(acf(logAQI[seq(1, n(), 7)], plot = F)$acf[2:4])),
    acf10.mean3 = mean(abs(acf(logAQI[seq(1, n(), 10)], plot = F)$acf[2:4])),
    .groups = "drop"
  )


### main =========================================================

# # read all data
dat <- do.call(
  rbind, lapply(yearrange, \(y) {
    f <- file.path(datapath, paste0(filenameprefix, y, ".csv"))
    print(f)
    dat <- read.csv(f)
    dat <- dat |>
      crop_latlon() |>
      left_join(
        siteinfo |>
          dplyr::select(
            Latitude, Longitude,
            Latitude.Binned, Longitude.Binned
          ),
        by = c("Latitude", "Longitude")
      ) |>
      dplyr::select(Latitude.Binned, Longitude.Binned, Date.Local, AQI) |>
      mutate(Date.Local = ymd(Date.Local)) |>
      filter(as.numeric(Date.Local - ymd(paste0(y, "-01-01"))) %% 7 == 0) |>
      group_by(Latitude.Binned, Longitude.Binned, Date.Local) |>
      summarise(AQI = mean(AQI), .groups = "drop") |>
      arrange(Date.Local)
    dat
  })
)

date_remove <- dat |>
  group_by(Date.Local) |>
  summarise(Nobs = n()) |>
  filter(Nobs == 1) |>
  pull(Date.Local)
dat <- dat |> filter(!(Date.Local %in% date_remove))

save(dat, file = file.path(datapath, "aqi-us.Rda"))


ggplot() +
  geom_sf(data = us_states) +
  geom_point(
    aes(
      x = Longitude, y = Latitude, color = Nrecord,
      size = Nrecord, group = Nrecord
    ),
    stroke = 1,
    data = dat |> filter(Date.Local == tail(Date.Local, 1)), shape = 1
  ) +
  theme_bw() +
  scale_color_continuous(type = "viridis") +
  labs(color = "Num Records", size = "Num Records")



load(file.path(datapath, "aqi-us.Rda"))



## visualize some subjects ===============
sampledat <- do.call(
  rbind, lapply(2022, \(y) {
    f <- file.path(datapath, paste0(filenameprefix, y, ".csv"))
    print(f)
    dat <- read.csv(f)
    dat <- dat |>
      crop_latlon() |>
      left_join(
        siteinfo |>
          dplyr::select(
            Latitude, Longitude,
            Latitude.Binned, Longitude.Binned
          ),
        by = c("Latitude", "Longitude")
      ) |>
      dplyr::select(
        Latitude, Longitude, Latitude.Binned, Longitude.Binned,
        Date.Local, AQI
      ) |>
      mutate(Date.Local = ymd(Date.Local)) |>
      arrange(Date.Local)
    dat
  })
)

selected_dates <- unique(sampledat$Date.Local)[c(1, 304)]

plots <- lapply(
  selected_dates, \(d) {
    ggplot() +
      geom_sf(data = us_states) +
      geom_point(aes(x = Longitude, y = Latitude),
        stroke = 1,
        data = sampledat |> filter(Date.Local == d),
        shape = 1, color = "#0B4590", alpha = 0.5
      ) +
      theme_bw() +
      labs(x = NULL, y = NULL, title = d) +
      theme(plot.title = element_text(hjust = 0.5))
  }
)

fig <- cowplot::plot_grid(plotlist = plots)

ggsave(
  filename = file.path("application", "api-sample.pdf"),
  fig, width = 8, height = 3,
  device = "pdf"
)


dat <- dat |>
  left_join(
    siteinfo |>
      dplyr::select(Latitude, Longitude, Latitude.Binned, Longitude.Binned),
    by = c("Latitude", "Longitude")
  )


sites <- dat |> distinct(Latitude.Binned, Longitude.Binned)

# rescaling latitude and longitude
sites$s <- (sites$Latitude.Binned - latrange[1]) / diff(latrange)
sites$t <- (sites$Longitude.Binned - lonrange[1]) / diff(lonrange)
dat$s <- (dat$Latitude.Binned - latrange[1]) / diff(latrange)
dat$t <- (dat$Longitude.Binned - lonrange[1]) / diff(lonrange)

# log transform of AQI
dat$y <- log1p(dat$AQI)
# par(mfrow=c(1,2))
# hist(dat$AQI, xlab="AQI", main="", breaks = 30)
# hist(dat$y, xlab="log(AQI+1)", main="", breaks = 30)
# par(mfrow=c(1,1))


## Settings ==============================
latgrid <- seq(latrange[1], latrange[2], 2)
longrid <- seq(lonrange[1], lonrange[2], 2)
locGrid <- margins2grid(list(latgrid, longrid))

evalGridList <- list(
  (latgrid - latrange[1]) / diff(latrange),
  (longrid - lonrange[1]) / diff(lonrange)
)
evalGrid <- margins2grid(evalGridList)

basis <- TensorBasis(list(
  create.bspline.basis(c(0, 1), nbasis = 6, norder = 4),
  create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
))
p <- attr(basis, "nbasis")

# B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)


## mean estimation ====================
fit.mean <- mgcv::bam(y ~ te(s, t), data = dat)
sites$y.mean <- mu <- as.numeric(predict(fit.mean, newdata = sites))

ggplot() +
  geom_sf(data = us_states) +
  geom_point(
    aes(
      x = Longitude.Binned, y = Latitude.Binned, color = y.mean,
      size = y.mean
    ),
    stroke = 1,
    data = sites, shape = 1
  ) +
  theme_bw() +
  scale_color_continuous(type = "viridis") +
  labs(
    x = "Longitude", y = "Latitude",
    color = "Average log(AQI+1)", size = "Average log(AQI+1)"
  )

dat$z <- dat$y - as.numeric(predict(fit.mean))

### Initialization ====================

# cov est function
mOpCov_path <- "external_codes/mOpCov/"
source(paste0(mOpCov_path, "mOpCov_prep.R"))
Rcpp::sourceCpp(paste0(mOpCov_path, "mOpCov_cpp.cpp"))

# read a sampled dataset
Lmi <- dat |>
  group_by(Date.Local) |>
  summarise(Nobs = n()) |>
  pull(Nobs)
dat$subjId <- rep(1:length(Lmi), Lmi)

sample_ratio <- 0.1
Lmi.sub <- if_else(Lmi < 10, Lmi, 10)

Ninit <- 100
set.seed(2501)
subjId.sub <- sort(sample(1:length(Lmi), Ninit, replace = FALSE))

dat.sub <- dat |> filter(subjId %in% subjId.sub)
ind <- lapply(subjId.sub, \(i) {
  ind <- rep(F, Lmi[i])
  idx <- sample(1:Lmi[i], Lmi.sub[i], replace = FALSE)
  ind[idx] <- TRUE
  return(ind)
}) |> unlist()
dat.sub <- dat.sub[ind, ]


locations <- dat.sub |>
  dplyr::select(s, t) |>
  as.matrix()
responses <- dat.sub$z
subject_ids <- rep(1:Ninit, Lmi.sub[subjId.sub])
q <- 4

CovRes <- mOpCov(
  location = locations, x = responses, subject = subject_ids,
  q = c(6, 6), lam = list(lam = 1e-10, alpha = 1e-6), ker = "cos"
)
FpcaOut <- fpca.mOpCov(OUT = CovRes)
lambdaInit <- FpcaOut$Eigen$values[1:q]
PhiInitEvalFull <- computeEigen(evalGrid, CovRes, FpcaOut)
PhiInitEval <- PhiInitEvalFull[, 1:q]

# # k <- 3
# ggplot() +
#   geom_sf(data=us_states) +
#   geom_point(aes(x=Longitude, y=Latitude, color=LogAQI, size=LogAQI), stroke=1,
#              data=data.frame(
#                Latitude=locGrid[,1], Longitude=locGrid[,2],
#                LogAQI=PhiInitEval[,k]
#              ), shape=1) +
#   theme_bw() +
#   ggtitle(paste("Eigenfunction", k)) +
#   scale_color_continuous(type = "viridis") +
#   labs(color = "LogAQI", size = "LogAQI")


ThetaInit <- smooth_basis(evalGrid, PhiInitEval, basis, lambda = 1e-8)$coefs
norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- eval_fd(evalGrid, FuncData(ThetaInit, basis))

# # fit diagonal element, estimate sigma2
# LtIdInit <- unlist(dat$Ltid[1:Ninit])
# innerId <- LtIdInit >= 0.25*max(LtIdInit) & LtIdInit <= 0.75*max(LtIdInit)
# covInitEvalDiag <- as.vector(PhiInitEvalFull^2 %*% FpcaOut$Eigen$values)
# sigma2Init <- mean((Yinit^2 - covInitEvalDiag[(LtIdInit)])[innerId])
# if (sigma2Init <= 1e-3) sigma2Init <- 1e-3
