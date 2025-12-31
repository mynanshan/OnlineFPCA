## ▶️ Application: files, outputs, and real data formats

All following scripts and results are stored in `application/`.

*NOTE*: By default, we set the root directory as the working directory instead of `application/`

### Essential application scripts

- `gfr.R` – runs the GFR data analysis pipeline (contains SBATCH headers and helper options). Produces `result_gfr*.Rdata` objects and supporting output for plotting and diagnostics.
- `analyze_gfr.R` – loads `result_gfr*.Rdata` and generates final figures and LaTeX tables (PDF :w
/ EPS / RDS outputs).
- `aqi.R` – runs the AQI data pipeline (2D spatial FPCA; contains SBATCH headers). Produces `fit_aqi*.Rdata` and batch comparisons. The commented code helps to produce `aqi_init.Rdata` that stores initializations.
- `analyze_aqi.R` – loads `fit_aqi*.Rdata` and generates maps / PDFs / RDS outputs for eigenfunctions and eigenvalues.
- `aqi_init.Rdata` (created once) – initial estimations used by AQI scripts.
- The script `aqi-eda.R` provides a complete walk-through of the exploratory data analysis of the AQI data and produces the plot `aqi-sample.pdf`.

### Result files (what to expect and their formats)

- `result_gfr.Rdata` – core OnlineFPCA output for the GFR application (contains objects like `PhiAvgEval`, `PhiAvgAll`, `lambda.avg`, `sgd_time`, `tau.select`, `tau_path`).
- `result_gfr_uq.Rdata` – OnlineFPCA output with uncertainty quantification (e.g. `resCI` / `CI` objects for CIs).
- `result_gfr_olcov.Rdata` – results from the competing OnlineCov method (eigenfunctions, timings)
- `fit_aqi.Rdata` – main AQI fit (OnlineFPCA results for AQI grid)
- `fit_aqi_uq.Rdata` – AQI fit with UQ (confidence intervals / CI objects)
- `fit_aqi_batch.Rdata` – batch/soap baseline results for AQI
- Figures – PDF/EPS files such as:
  - `gfr_eigfun.pdf`, `gfr_taupath.pdf`
  - `aqi-eigfun.pdf`, `aqi-taupath.pdf` etc.
- RDS artifacts (saved for later use):
  - `eigf_gfr.rds`, `eigval_gfr.rds`
  - `eigf_aqi.rds`, `eigval_aqi.rds`, `tgrid_aqi.rds`

Each `.Rdata` usually contains named R objects used by `analyze_*` scripts; `.rds` files store single R objects (e.g. evaluated eigenfunctions or eigenvalues), and PDFs/EPS are directly used in the manuscript figures.

---

## ▶️ Reproducing application results (real data)

Each application has two main steps:
1. Run the heavy computation pipeline that estimates models and saves `*.Rdata` objects to `` (may run on a cluster).
2. Run the corresponding `analyze_*` script that loads the `*.Rdata` files and produces figures/tables (PDF / EPS / LaTeX tables).

Examples:

```bash
# Run the GFR application (single run or on cluster, this script contains SBATCH headers)
Rscript application/gfr.R

# After the run finishes, generate the publication-ready figures
Rscript application/analyze_gfr.R
# output: files like `gfr_eigfun.pdf`, `gfr_taupath.pdf`, `gfr_time.eps`

# Similarly for AQI
Rscript application/aqi.R
Rscript application/analyze_aqi.R
# outputs saved under `` (e.g. `aqi_<figure>.pdf`)
```

---

### Real data sources & expected formats

- GFR data (`data/gfr/`):
  - **`GFR2023Feb01.csv`** – the main tabular dataset used by `gfr.R`. The script expects per-subject repeated measurements in (the default dataset) columns 3:9 (one column per year) and uses these to build `Ly` and `Ltid` lists.
  - The script performs simple preprocessing (removes subjects with all-NA records, rescales times, de-mean using column means).
  - The original data is only publicly available. We have attached a `GFR-pseudo.csv` as a pseudo dataset for illustrative purposes. 

- AQI data (`data/epa-aqs/`):
  - **`aqi-us.Rda`** – a preprocessed dataset used directly by `aqi.R` and `analyze_aqi.R`.
  - **Daily CSVs** (e.g. `daily_81102_1982.csv`, ... `daily_81102_2022.csv`) – the raw source files used to create the summarized `aqi-us.Rda`.
  - **`siteinfo.xlsx`** – site metadata (locations) used to create `locGrid` and mapping visualizations.
  - Expected fields in the AQI data include: `Latitude.Binned`, `Longitude.Binned`, `Date.Local`, `AQI` (used to compute `y = log1p(AQI)`).
