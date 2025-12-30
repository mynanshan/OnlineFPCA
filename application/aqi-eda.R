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

# Preprocessing script for AQI data
# - Builds `siteinfo.xlsx` (site-level metadata with binned coordinates)
# - Reads raw daily AQI CSVs, crops and bins locations, subsamples weekly,
#   aggregates by (binned location, date), removes singleton-date observations,
#   and saves `data/epa-aqs/aqi-us.Rda`.
# - Exports a sample figure `application/api-sample.pdf` illustrating the data.
#
# Usage: `Rscript application/aqi-eda.R`
# Required packages (non-exhaustive): readr, dplyr, tidyr, lubridate, giscoR, sf,
# spData, ggplot2, openxlsx

# source("./R/fdaMdim.R")
# source("./R/helper.R")

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


# Helper functions
# crop_latlon: keep only points inside the configured lon/lat bounding box
crop_latlon <- function(dat) {
  dat |>
    filter(Longitude > lonrange[1] &
      Longitude < lonrange[2] &
      Latitude > latrange[1] &
      Latitude < latrange[2])
}

# bin_latlon: bin continuous lat/lon values to a regular grid of width `width`.
# Returns additional columns `Latitude.Binned` and `Longitude.Binned` by default.
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

# Quick interactive checks removed (kept minimal, reproducible plots below)

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


## Diagnostics (optional)
# The following diagnostics (ACF and small checks) are useful during
# exploration but are disabled by default to keep this preprocessing script
# deterministic and fast. Set `if (TRUE)` to run them interactively.
if (FALSE) {
  # Check ACF on 2022 data (subsampled and binned)
  dat <- read.csv(file.path(datapath, "daily_81102_2022.csv"))
  dat <- dat |>
    crop_latlon() |>
    bin_latlon() |>
    select(Latitude.Binned, Longitude.Binned, Date.Local, AQI) |>
    group_by(Latitude.Binned, Longitude.Binned, Date.Local) |>
    summarise(AQI = mean(AQI), .groups = "drop") |>
    mutate(Date.Local = ymd(Date.Local)) |>
    arrange(Date.Local)
  # TODO: consider geometric mean (mean of log(AQI)) as alternative aggregation.

  dat |>
    group_by(Date.Local) |>
    summarise(Nobs = n())

  # Example plot for a single date
  ggplot() +
    geom_sf(data = us_states) +
    geom_point(aes(x = Longitude.Binned, y = Latitude.Binned, color = log1p(AQI)),
      stroke = 1,
      data = dat |> filter(Date.Local == ymd("2022-03-01")), shape = 1
    ) +
    theme_bw() +
    scale_color_continuous(type = "viridis") +
    labs(color = "AQI")

  # Example ACF computations
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
}


### main =========================================================
# Read yearly raw CSVs, crop and bin by location, subsample weekly (one observation
# per 7-day stride from Jan 1), aggregate duplicated grid/day pairs by mean AQI,
# and remove dates with only a single observation across all sites.

# Read and aggregate per-year files
dat <- do.call(
  rbind, lapply(yearrange, \(y) {
    f <- file.path(datapath, paste0(filenameprefix, y, ".csv"))
    message("Reading ", f)
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
      # Subsample weekly anchored at Jan 1 of that year
      filter(as.numeric(Date.Local - ymd(paste0(y, "-01-01"))) %% 7 == 0) |>
      # Aggregate in case multiple monitors fall in the same bin/date
      group_by(Latitude.Binned, Longitude.Binned, Date.Local) |>
      summarise(AQI = mean(AQI), .groups = "drop") |>
      arrange(Date.Local)
    dat
  })
)

# Remove dates that have only one observation across all locations (not informative)
date_remove <- dat |>
  group_by(Date.Local) |>
  summarise(Nobs = n()) |>
  filter(Nobs == 1) |>
  pull(Date.Local)
dat <- dat |> filter(!(Date.Local %in% date_remove))

# Save processed data for downstream analyses
save(dat, file = file.path(datapath, "aqi-us.Rda"))
message("Saved ", file.path(datapath, "aqi-us.Rda"), " (rows=", nrow(dat), ", dates=",
  paste(range(dat$Date.Local), collapse = " to "), ")")

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
