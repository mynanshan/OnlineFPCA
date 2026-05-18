## Application: Files, Outputs, and Real Data Formats

All following scripts and results are stored in `application/`.

*NOTE*: By default, we set the root directory as the working directory instead of `application/`.

### Essential application scripts

- `application/gfr/gfr.R` - runs the GFR data analysis pipeline. Produces `result_gfr*.Rdata` objects and supporting output for plotting and diagnostics.
- `application/gfr/analyze_gfr.R` - loads `result_gfr*.Rdata` and generates final figures and LaTeX tables (PDF / EPS / RDS outputs).
- `application/aqi/aqi.R` - runs the AQI data pipeline (2D spatial FPCA). Produces `fit_aqi*.Rdata` and batch comparisons. The commented code helps produce `aqi_init.Rdata`, which stores initializations.
- `application/aqi/analyze_aqi.R` - loads `fit_aqi*.Rdata` and generates maps, PDFs, and RDS outputs for eigenfunctions and eigenvalues.
- `application/aqi/aqi-eda.R` - provides a complete walk-through of the exploratory data analysis of the AQI data and produces `aqi-sample.pdf`.
- `application/news/news-eda.R` - preprocesses the online news data from `data/news/Data/` into `data/news/processed/news_preprocessed.rds` and `news_preprocess_summary.csv`.
- `application/news/news.R` - runs the news-click application. It supports `--run fpca`, `--run ci`, `--run olcov`, or `--run all`, and produces `result_news*.Rdata` objects.
- `application/news/analyze_news.R` - loads the news fits and preprocessing object, then generates final news figures, summaries, and RDS artifacts.
- `aqi_init.Rdata` (created once) - initial estimations used by AQI scripts.

### Result files

- GFR:
  - `result_gfr.Rdata` - core OnlineFPCA output for the GFR application, including objects such as `PhiAvgEval`, `PhiAvgAll`, `lambda.avg`, `sgd_time`, `tau.select`, and `tau_path`.
  - `result_gfr_uq.Rdata` - OnlineFPCA output with uncertainty quantification, such as `resCI` / `CI` objects.
  - `result_gfr_olcov.Rdata` - results from the competing OnlineCov method.
  - Figures include `gfr_eigfun.pdf`, `gfr_taupath.pdf`, and timing outputs such as `gfr_time.eps`.
  - RDS artifacts include `eigf_gfr.rds` and `eigval_gfr.rds`.
- AQI:
  - `fit_aqi.Rdata` - main AQI fit (OnlineFPCA results for AQI grid).
  - `fit_aqi_uq.Rdata` - AQI fit with UQ (confidence interval objects).
  - `fit_aqi_batch.Rdata` - batch / SOAP baseline results for AQI.
  - Figures include `aqi-eigfun.pdf`, `aqi-taupath.pdf`, and related map outputs.
  - RDS artifacts include `eigf_aqi.rds`, `eigval_aqi.rds`, and `tgrid_aqi.rds`.
- News:
  - `result_news.Rdata` - main OnlineFPCA fit for the news-click application.
  - `result_news_uq.Rdata` - OnlineFPCA fit with uncertainty quantification.
  - `result_news_olcov.Rdata` - OnlineCov comparison fit.
  - Figures and summaries include `news_eigfun.pdf`, `news_mean.pdf`, `news_taupath.pdf`, `news_time.pdf`, `news_fve.csv`, and `news_summary.txt`.
  - RDS artifacts include `eigf_news.rds` and `eigval_news.rds`.

Each `.Rdata` usually contains named R objects used by `analyze_*` scripts. Each `.rds` file stores a single object, such as evaluated eigenfunctions or eigenvalues. PDFs, EPS files, CSV summaries, and text summaries are directly used for manuscript figures and reporting.

---

## Reproducing Application Results

Each application has two main steps:

1. Run the heavy computation pipeline that estimates models and saves `*.Rdata` objects under `application/` (often on a cluster).
2. Run the corresponding `analyze_*` script that loads the `*.Rdata` files and produces figures, tables, and reusable RDS artifacts.

Examples:

```bash
# GFR
Rscript application/gfr/gfr.R
Rscript application/gfr/analyze_gfr.R

# AQI
Rscript application/aqi/aqi.R
Rscript application/aqi/analyze_aqi.R

# News preprocessing, fitting, and analysis
Rscript application/news/news-eda.R
Rscript application/news/news.R --run all
Rscript application/news/analyze_news.R
```

---

## Real Data Sources and Expected Formats

- GFR data (`data/gfr/`):
  - `GFR2023Feb01.csv` - the main tabular dataset used by `application/gfr/gfr.R`. The script expects per-subject repeated measurements in columns 3:9 (one column per year) and uses these to build `Ly` and `Ltid` lists.
  - The script performs simple preprocessing: removes subjects with all-NA records, rescales times, and de-means using column means.
  - The original data is not publicly available. The repository includes `GFR-pseudo.csv` as a pseudo dataset for illustrative purposes.

- AQI data (`data/epa-aqs/`):
  - `aqi-us.Rda` - a preprocessed dataset used directly by `application/aqi/aqi.R` and `application/aqi/analyze_aqi.R`.
  - Daily CSVs such as `daily_81102_1982.csv`, ..., `daily_81102_2022.csv` - raw source files used to create the summarized `aqi-us.Rda`.
  - `siteinfo.xlsx` - site metadata used to create `locGrid` and mapping visualizations.
  - Expected fields include `Latitude.Binned`, `Longitude.Binned`, `Date.Local`, and `AQI`; the application computes `y = log1p(AQI)`.

- News data (`data/news/`):
  - `data/news/Data/News_Final.csv` - metadata for news items.
  - `data/news/Data/<Platform>_<Topic>.csv` - social feedback trajectories for Facebook, GooglePlus, and LinkedIn topics.
  - `data/news/processed/news_preprocessed.rds` - preprocessed functional data object produced by `application/news/news-eda.R`.
  - `data/news/processed/news_preprocess_summary.csv` - preprocessing counts and retention summary produced by `application/news/news-eda.R`.
