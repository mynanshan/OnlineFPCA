# OnlineFPCA
A compact, reproducible codebase for the OnlineFPCA project used in the associated paper.

This repository contains R code to run online stochastic-gradient FPCA, baseline batch
methods, simulation studies, and real-data analyses for the paper. The README below
explains the project structure, how to set up the environment, and how to reproduce
tables and figures from the paper.  

**NOTE**: the `output/` directory contains a large number of experiment result files, which may slow down repo cloning.

---

## 📦 Repository structure

Top-level layout (important files/directories):

- **`R/`** – Core R implementation files
    - `onlineFPCA.R`: OnlineFPCA algorithm
    - `onlineFDAlocalpoly.R`: OnlineCov algorithm (adapted from the code for Yang & Yao (2023))
    - `fpcaReg.R`: SOAP algorithm, a batch FPCA method
    - Other helper utilities.
- **`experiments/`** – Simulation scripts and analysis helpers (e.g. `fpca1d.R`, `fpca2d.R`). `run_reps.sh` is a helper for submitting repeated (seeded) simulation jobs on an HPC cluster.
- **`test/`** – `demo-fpca1d.R` and `demo-fpca2d.R` are also for 1D and 2D simulation studies, but also easier to be run interactively. 
- **`application/`** – Real-data pipelines and plotting scripts (e.g. `gfr.R`, `aqi.R`).
- **`data/`** – Data required to reproduce the application results. For this project the subfolders of interest are `data/gfr/` and `data/epa-aqs/`.
- **`data_generation/`** – Data simulators used by the simulation studies.
- **`external_codes/`** – Third-party algorithms: mOpCov and REML.
- **`install_pkgs.R`** – Script to install required R packages in a fresh environment.
- **`output/`** – Include our result files from simulation studies and real-data analysis. Large files are removed due to size limits.

---


## 🧪 Reproduce tables and figures

**Main article**:

- Table 1, 2, S4: Main simulation studies for 1D and 2D data under different noise levels. In `experiments/`, Run `fpca1d.R` and `fpca2d.R` for `--seed` across 1~100. Then, run `analyze_fpca1d.R` and `analyze_fpca2d.R` to obtain results from `knitr::kable()`. Table 1 and 2 are part of Table S4.

- Figure 3 and S6: Dynamic tuning path in 1D/2D simulation. The script `fpca1d.R` will produce the dynamic tuning path `taupath-sim*d-sgd.pdf` when the seed is `1234`. You can set `simple` to `1` and `compare` to `0` to simplify the workflow and focus on figure plotting. 
```r
Rscript experiments/fpca1d.R --seed 1234 --compare 0 --simple 1
Rscript experiments/fpca2d.R --seed 1234 --compare 0 --simple 1
```

- Figure 4: Subject examples from the AQI data. The script `application/aqi-eda.R` provides a complete walk-through of the exploratory data analysis of the AQI data, and the command to generate Figure 4 (`aqi-sample.pdf`). 

- Figure 5: FPC plot for AQI data. First run `application/aqi.R` to obtain `fit_aqi_*.Rdata`. Then, run `application/analyze_aqi.R` to produce the plot `simu-aqi-fpc.pdf`.

- Figure 6: FPC plot for GFR data. First run `application/gfr.R` to obtain `result_gfr_*.Rdata`. Then, run `application/analyze_gfr.R` to produce the plot `gfr_eigfun.pdf`.

- Figure S1: Use `Matrix::image()` to plot the matrix `B`, `G` and `S` in the 1D simulation.

- Figure S2: CI plot for 1D simulation. In `experiments/`, run `ci-simu1d.R` for `--seed` across 1~100. Then, run `analyze_ci1d.R` to produce the plot `ci1d-ci.pdf`. 

- Figure S3: CI plot for 2D simulation. In `experiments/`, run `ci-simu2d.R` for `--seed` across 1~100. Then, run `analyze_ci2d.R` to produce the plot `ci2d-ci.pdf`. 

- Table S1: RMSEs for different $(C,W,B)$. Run `check_bwc.R` for `--seed` across 1~100. Then, `analyze_bwc.R` will produce the table.

- Figure S4, Table S2: RMSEs for different $\omega$ and block sizes. Run `check_abv.R` for `--seed` across 1~100. Then, run `analyze_abv.R` to produce the figure `abv_paths.pdf` and the table.

- Table S3: RMSEs for different mini-batch sizes.  Run `check_nbatch.R` for `--seed` across 1~100. Then, run `analyze_nbatch.R` to produce the table.


---

## 🧰 Environment setup

1. Recommended R (tested): R 4.5.0 (or a recent patch release).

2. Install the packages used in the analyses:

```bash
Rscript install_pkgs.R
```

3. Make sure `data/` contains the required inputs. The (large) application datasets are not stored in the repo; if needed, place the original data files in `data/gfr/` and `data/epa-aqs/` as used by `application/gfr.R` and `application/aqi.R`.

---

## Tested Package Versions

| Package | Version |
|---|---|
| argparse | 2.3.1 |
| cowplot | 1.2.0 |
| doFuture | 1.1.3 |
| dplyr | 1.1.4 |
| face | 0.1-8 |
| fastmatrix | 0.6-4 |
| fda | 6.3.0 |
| fdapace | 0.6.0 |
| foreach | 1.5.2 |
| giscoR | 1.0.0 |
| gslnls | 1.4.2 |
| lubridate | 1.9.4 |
| ManifoldOptim | 1.0.1 |
| Metrics | 0.1.4 |
| mgcv | 1.9-4 |
| Rcpp | 1.1.0 |
| Rdimtools | 1.1.3 |
| readr | 2.1.6 |
| readxl | 1.4.5 |
| remotes | 2.5.0 |
| RSpectra | 0.16-2 |
| rTensor | 1.4.9 |
| sf | 1.0-23 |
| sm | 2.2-6.0 |
| spData | 2.3.4 |
| stringr | 1.6.0 |
| tidyr | 1.3.2 |