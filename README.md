# OnlineFPCA
A compact, reproducible codebase for the OnlineFPCA project used in the associated paper.

This repository contains R code to run online stochastic-gradient FPCA, baseline batch
methods, simulation studies, and real-data analyses for the paper. The README below
explains the project structure, how to set up the environment, and how to reproduce
tables and figures from the paper.  

---

## 📦 Repository structure

Top-level layout (important files/directories):

- **`R/`** – Core R implementation files
    - `onlineFPCA.R`: OnlineFPCA algorithm
    - `onlineFDAlocalpoly.R`: OnlineCov algorithm (adapted from the code for Yang & Yao (2023))
    - `fpcaReg.R`: SOAP algorithm, a batch FPCA method
    - Other helper utilities.
- **`experiments/`** – Simulation scripts and analysis helpers (e.g. `fpca1d.R`, `fpca2d.R`).
- **`application/`** – Real-data pipelines and plotting scripts (e.g. `gfr.R`, `aqi.R`).
- **`data/`** – Data required to reproduce the application results. For this project the subfolders of interest are `data/gfr/` and `data/epa-aqs/`.
- **`data_generation/`** – Data simulators used by the simulation studies.
- **`external_codes/`** – Third-party algorithms: mOpCov and REML.
- **`install_pkgs.R`** – Script to install required R packages in a fresh environment.
- **`experiments/run_reps.sh`** – Helper for submitting repeated (seeded) simulation jobs on an HPC cluster.

---

## 🧰 Environment setup

1. Recommended R (tested): R 4.5.0 (or a recent patch release).

2. Install the packages used in the analyses:

```bash
Rscript install_pkgs.R
```

3. Make sure `data/` contains the required inputs. The (large) application datasets are not stored in the repo; if needed, place the original data files in `data/gfr/` and `data/epa-aqs/` as used by `application/gfr.R` and `application/aqi.R`.

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