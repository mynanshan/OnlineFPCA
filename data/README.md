# Data

## EPA Daily AQI Data

**Data download.** [[Link](aqs.epa.gov/aqsweb/airdata/download_files.html)] Select the files in the daily PM10 data (series 81102). Downloaded files are named `daily_81102_yyyy.csv`, where `yyyy` stands for the year. 

**Data document.** [[Link](https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_content_4)] Contains variable names and meaning in the raw data files.

**Preprocessing.**
The preprocessing steps below are implemented in `application/aqi-eda.R` to produce `aqi-us.Rda` that aggregates all PM10 data and `siteinfo.xlsx` that stores all observation locations. To run the application study, one simply need to load the preprocessed data. Below summarizes the steps performed by `application/aqi-eda.R`:

1. Spatial / temporal bounds and grid
   - Longitude range: **-125 to -67**, Latitude range: **25 to 49**.
   - Grid width: **2°** (cells centered at half-grid offsets). These define the binned location grid used in later analyses.

2. Build `siteinfo` (site-level metadata)
   - Read per-year raw files and aggregate records grouped by (State.Code, County.Code, Site.Num).
   - For each site, compute unique Latitude/Longitude and total `Nrecord` across years.
   - Crop sites to the bounding box and bin their lat/lon to the 2° grid.
   - For sites that do not fall within the U.S. polygon, project them to the nearest in-grid location inside the U.S. mask.
   - Save metadata to `data/epa-aqs/siteinfo.xlsx` (contains e.g. original Latitude/Longitude, binned coordinates, and `Nrecord`).

3. Aggregate daily observations into the grid and weekly sampling
   - For each year `yyyy`, read `daily_81102_yyyy.csv`, crop to bounding box and join with `siteinfo` to obtain binned coordinates.
   - Select and keep: `Latitude.Binned`, `Longitude.Binned`, `Date.Local`, and `AQI`.
   - Convert `Date.Local` to `Date` objects and subsample **weekly** by keeping dates where `(Date.Local - Jan1) %% 7 == 0` (i.e., one observation per 7-day stride from Jan 1 of each year).
   - For each (binned location, date) take the **mean AQI** (aggregate duplicates if multiple monitors fall into the same bin/date).

4. Clean up and save
   - Remove days that have only **one** observation across all locations (these dates provide no useful spatial comparison).
   - Save the processed dataset as `data/epa-aqs/aqi-us.Rda`. The saved object `dat` contains at least the columns:
     - `Latitude.Binned`, `Longitude.Binned` (binned grid coordinates)
     - `Date.Local` (as `Date`)
     - `AQI` (mean AQI for that grid cell and date)

5. Notes for downstream analyses
   - Downstream scripts (e.g., `application/aqi.R`) perform additional steps such as
     - computing `y = log1p(AQI)` (log-transform of AQI),
     - rescaling binned lat/lon to [0,1] (`s,t`) for the model grid,
     - building basis representations for spatial FPCA.
   - These transformations are applied at analysis time and are *not* required in the saved `aqi-us.Rda` file.


## GFR Kidney Transplant Data

**Access.** The GFR dataset used in this project is private and must be requested from the data custodian. See the following paper for data request instructions.
* Liu, H., You, J., & Cao, J. (2023). Functional L-optimality subsampling for functional generalized linear models with massive data. Journal of Machine Learning Research, 24(219), 1-41.


**Source file:** `data/gfr/GFR2023Feb01.csv` (not included in this repository. We have attached a `GFR-pseudo.csv` as a pseudo dataset for illustrative purposes.).

**File format and contents:**
- Rows: one row per subject (e.g. 144,703 rows in the provided CSV).
- Columns (9 total):
  1. **ID** — subject identifier (unique per patient).
  2. **Y** — binary indicator (1 if patient survives more than 10 years after transplant, 0 otherwise).
  3–9. **GFR measurements** — annual GFR values for years 1..7 after transplant (columns 3 to 9). Missing values are allowed (NA) when a patient lacks a measurement for a given year.

**Preprocessing used by our analysis (`application/gfr.R`):**
1. Remove subjects with no GFR observations (rows where all columns 3:9 are NA).
   - Implemented as: `rm_idx <- apply(is.na(dat[, 3:9]), 1, sum) >= 7; dat <- dat[!rm_idx, ]`.
2. De-mean measurement vectors by column means (computed with `na.rm = TRUE`):
   - `gfr_mean <- colMeans(as.matrix(dat[, 3:9]), na.rm = TRUE)`
   - For each subject i, the observed vector is `dat$Ly[[i]] = Yi - B[Ti, ] %*% theta_mu` in downstream code; in preprocessing we compute `dat$Ly` as the subject's observed GFR values with the global column means subtracted.
3. Build per-subject lists used by the model:
   - `Ly` — list of observed values for each subject (only non-NA entries),
   - `Ltid` — corresponding time indices (integer positions of observed years; these are then converted to a scaled time grid used in models),
   - `Lmi` — number of observations per subject.
   - In the code this is done by mapping non-missing entries of columns 3:9 into `Ly` / `Ltid` lists.
4. Rescale time indices to [0,1] (used by FPCA routines): `Lt = (ti - 1) / (nt - 1)` where `nt` is the total number of grid points used.
