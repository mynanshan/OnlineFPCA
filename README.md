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

- Table 1, 2: Main simulation studies for 1D and 2D data. In `experiments/`, Run `fpca1d.R` and `fpca2d.R` for `--seed` across 1~100. Then, run `analyze_fpca1d.R` and `analyze_fpca2d.R` to obtain results from `knitr::kable()`.

- Figure 3: Dynamic tuning path in 1D simulation. The script `fpca1d.R` will produce the dynamic tuning path when the seed is `1234`. You can set `simple` to `1` and `compare` to `0` to simplify the workflow and focus on figure plotting. 
```r
Rscript experiments/fpca1d.R --seed 1234 --compare 0 --simple 1
```

- Figure 4: Subject examples from the AQI data. The script `application/aqi-eda.R` provides a complete walk-through of the exploratory data analysis of the AQI data, and the command to generate Figure 4 (`aqi-sample.pdf`).  
