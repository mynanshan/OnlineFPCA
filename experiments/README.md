# Simulation Studies

All following scripts and results are stored in `experiments/`.

*NOTE*: By default, we set the root directory as the working directory instead of `experiments/`

Experiment files:

* `fpca1d.R`: main 1D simulation
* `fpca2d.R`: main 2D simulation
* `ci-simu1d.R`: 1D simulation with pointwise CI
* `ci-simu2d.R`: 2D simulation with pointwise CI
* `check_abv.R`, `plot_abv.R`: Check how $\omega$ and block size affects 1D simulation performance
* `check_bwc.R`: Check how $(C,W,B)$ affects 1D simulation performance
* `check_nbatch.R`: Check how batch size affects 1D simulation performance
* `size1d.R`: compare FPC RMSEs across sample sizes for 1D simulations
* `size2d.R`: compare FPC RMSEs across sample sizes for 2D simulations
* `simu_gfr.R`: simulation using FPCs from GFR data
* `simu_aqi.R`: simulation using FPCs from AQI data

Result analysis:

* `analyze_<exprmt-id>.R`

## Results (formats and what the analysis scripts produce)

This project follows a consistent output pattern for simulation runs and subsequent analyses. Result files are generally in a folder labelled as the experiment name. Below are the common file formats and examples for a few important experiments:

- Per-seed CSV outputs
  - Naming: `experiments/<exprmt>/simu_<exprmt>_sd<seedtext>.csv` (where `<seedtext>` is a zero-padded seed id).
  - Contents: a row per checkpoint (or iteration) and columns typically include metadata and performance metrics, for example:
    - `seed`, `Method`, `StepSize`, `N`, `Ninit`, `initMethod`, `npc`, `nBatch`, `Time`, `RMSEphi1`, `RMSEphi2`, `RMSEphi3`, `RMSEphi1.avg`, ...
  - Example: `fpca1d.R` and `simu_gfr.R` produce one CSV file per seed storing RMSEs and run-time information.

- Per-seed RDS/RData outputs (optional)
  - Naming: often `Theta_<exprmt>_sd<seedtext>.rds` or `ci_<exprmt>_sd<seedtext>.RData` for CI-related outputs.
  - Contents: richer objects used by analyze scripts (e.g., `ThetaRecord` lists of averaged `Theta` matrices, or `CIobj` containing confidence intervals).
  - Example: `simu_gfr.R` and `simu_aqi.R` save `Theta_<exprmt>_sd*.rds` (per-seed FPC records), and `ci-simu1d.R`/`ci-simu2d.R` save `CIobj` RData files.

- What the `analyze_<exprmt-id>.R` scripts do
  - Aggregate per-seed CSVs across seeds into a single table and compute summary statistics (means, standard errors, etc.).
  - Produce LaTeX tables (via `knitr::kable`) summarizing averaged metrics (e.g., average RMSEs across 100 repeats for `fpca1d`).
  - Read optional RDS/RData artifacts (e.g., `Theta_*.rds`, `ci_*.RData`) to compute and plot average eigenfunctions, confidence bands, and other diagnostics.
  - Save plots as PDF/EPS and summary tables to the experiment directory.
  - Example: `analyze_fpca1d.R` reads the per-seed CSVs and outputs a LaTeX table of average RMSEs (across all seeds); `analyze_simu_gfr.R` produces a summary table and FPC plots (`simu-gfr-fpc.pdf`).


## ▶️ Reproducing simulation studies

General workflow for simulation experiments:

1. Run per-seed simulation jobs that save per-seed CSV/RDS outputs to `experiments/<exprmt>/` (where `<exprmt>` is `fpca1d`, `fpca2d`, ...). Example (single seed):

```bash
# Run one seed locally
Rscript experiments/fpca1d.R --seed 1

# Or submit many seeds to an HPC scheduler using the helper (uses sbatch internally):
./experiments/run_reps.sh experiments/fpca1d.R 100
```

2. After the per-seed outputs are available, run the corresponding analysis/aggregate script to produce tables and figures. For example:

```bash
# Aggregate and generate tables/plots for the GFR simulation study
Rscript experiments/analyze_fpca1d.R
```

The analyze scripts read the `experiments/<exprmt>/` directory and produce tables/plots (see the top of `analyze_fpca1d.R` for the expected file names). Replace `<exprmt>` with the experiment name.

---

Below are **detailed**, per-experiment reproduction notes (what the script does, how to run it, and the expected output files). Use `--seed <int>` to run individual seeds; for large sweeps use `./experiments/run_reps.sh <script> <#seeds>` which calls the script with seeds 1..N.

---

### 1D Simulation ✅
- Purpose: point estimation of FPCs and comparison vs. batch methods; includes a CI variant for FPC CIs.
- Run (single seed / local):
  - Point estimation: `Rscript experiments/fpca1d.R --seed 1`
  - Aggregation/analysis: `Rscript experiments/analyze_fpca1d.R`
  - Plot a single tuning path (example): `Rscript experiments/fpca1d.R --seed 1234 --compare 0 --simple 1` (produces `taupath-sim1d-sgd.pdf` for the canonical seed)
- CI variant:
  - Collect per-seed CI runs: `./experiments/run_reps.sh experiments/ci-simu1d.R 100`
  - Aggregate/plot CIs: `Rscript experiments/analyze_ci1d.R` (produces `ci1d-ci.pdf`)
- Typical outputs (per seed):
  - `experiments/fpca1d/simu_fpca1d_sdXXXX.csv` — CSV rows per checkpoint with columns like `seed, Method, StepSize, N, Ninit, npc, nBatch, Time, RMSEphi1, RMSEphi2, ...`
  - `experiments/ci1d/ci_ci1d_sdXXXX.RData` — RData with `CIobj` (contains `PhiStar`, `PhiTrue`, and `CIs` arrays)
- Notes: use `--seed` for debugging one seed; CI runs save intermediate CI arrays and a global optimum computed via ManifoldOptim.

---

### 2D Simulation ✅
- Purpose: extend 1D experiments to 2D surfaces (same pipelines but tensor bases / mOpCov initialization).
- Run (single seed / local):
  - Point estimation: `Rscript experiments/fpca2d.R --seed 1` (flags: `--compare` to run batch baselines, `--simple` for quick runs)
  - Aggregation/analysis: `Rscript experiments/analyze_fpca2d.R`
- CI variant:
  - Per-seed CI runs: `./experiments/run_reps.sh experiments/ci-simu2d.R 100`
  - Aggregate/plot CIs: `Rscript experiments/analyze_ci2d.R` (produces `ci2d-ci.pdf`)
- Typical outputs (per seed):
  - `experiments/fpca2d/simu_fpca2d_sdXXXX.csv`
  - Example plot: `experiments/fpca2d/taupath-sim2d-<sgdtype>.pdf` (produced for canonical seeds)
  - `experiments/ci2d/ci_ci2d_sdXXXX.RData` — similar `CIobj` format as 1D
- Notes: 2D experiments rely on `external_codes/mOpCov/` for initialization and are more expensive than 1D.

---

### Confidence Interval Experiments (CI) 🧪
- Scripts: `ci-simu1d.R` and `ci-simu2d.R`
- Run: `Rscript experiments/ci-simu1d.R --seed <s>` or run via `run_reps.sh` for many seeds
- Output: `experiments/ci1d/ci_ci1d_sdXXXX.RData` (and `ci2d` analog), with `CIobj` containing fields:
  - `PhiStar`: numerically optimized FPCs
  - `PhiTrue`: the truth used in simulation
  - `CIs`: array of bootstrap/CI draws per FPC
- Analysis: `Rscript experiments/analyze_ci1d.R` / `analyze_ci2d.R` → CI plots (`ci1d-ci.pdf`, `ci2d-ci.pdf`) and numerical summaries
- Notes: CI runs set `fpcCI = TRUE` inside `fpca.sgd()` and compute global optima with ManifoldOptim for reference.

---

### Check mini-batch size (nbatch) ⚖️
- Script: `check_nbatch.R` (produces results for a grid of `nBatch` values)
- Run: `./experiments/run_reps.sh experiments/check_nbatch.R 100`
- Current output per seed: `experiments/abv1d/simu_abv1d_sdXXXX.csv` with columns like `seed, Method, epoch, nBatch, nBlock, tau, RMSEphi1, ...`
- Analysis: `Rscript experiments/analyze_nbatch.R` aggregates RMSEs by `nBatch` and produces LaTeX tables. Note that this analysis script currently reads from `experiments/nbatch1d/simu_nbatch1d_sdXXXX.csv`, so either generate/copy results into that layout or adjust the script before running the aggregation.
- Notes: This uses parallel execution (`doFuture`) to sweep over `nBatch` values faster.

---

### Sample-size experiments (size1d and size2d) 📈
- Purpose: compare FPC estimation RMSEs across increasing sample sizes.
- 1D run:
  - Single seed and method: `Rscript experiments/size1d.R --seed 1 --algo 1`
  - `--algo` choices: `1` OnlineFPCA, `2` PACE, `3` FACE, `4` REML, `5` SOAP
  - Repeat across seeds and algorithms before analysis.
- 2D run:
  - Single seed and method: `Rscript experiments/size2d.R --seed 1 --algo 1`
  - `--algo` choices: `1` OnlineFPCA, `2` SOAP, `3` mOpCov
  - Optional controls include `--N`, `--ninit`, `--nsmall`, `--mopcov-max-n`, and `--mopcov-max-obs`.
- Typical outputs:
  - `experiments/size1d/simu_size1d_alg*_sdXXXX.csv`
  - `experiments/size2d/simu_size2d_alg*_sdXXXX.csv`
- Analysis:
  - `Rscript experiments/analyze_size1d.R` → `experiments/size1d/size1d_rmse.pdf`
  - `Rscript experiments/analyze_size2d.R` → `experiments/size2d/size2d_rmse.pdf`

---

### Check dynamic tuning (ABV) and ABV plot 🔧
- Script: `check_abv.R` (sweeps `nBlock` and `ewmabv.beta` values)
- Run: `./experiments/run_reps.sh experiments/check_abv.R 100`
- Outputs:
  - `experiments/abv1d/simu_abv1d_sdXXXX.csv` — per-seed ABV results (epoch-wise RMSEs and chosen `tau`)
  - `experiments/plot_abv.R` can produce diagnostic ABV path figures (e.g., `abv_paths.pdf`)
- Analysis: `Rscript experiments/analyze_abv.R` (produces `abv_paths.pdf` and a summary table of RMSEs)
- Notes: ABV experiments are useful to tune the ABV exponential-weight parameter and block-size choices.

---

### Check BWC (C, W, tau selection grid) 🧭
- Script: `check_bwc.R` (iterates over a grid of `(C, W)` tuning parameters)
- Run: `./experiments/run_reps.sh experiments/check_bwc.R 100`
- Outputs:
  - `experiments/bwc/simu_bwc_sdXXXX.csv` — contains columns `C`, `W`, and the usual RMSE metrics
  - `experiments/bwc/bwcSet.csv` (saved when `seed==1`) documents the `(C,W)` grid used
  - Per-seed `taupath_<basename>.csv` entries saved for debugging
- Analysis: `Rscript experiments/analyze_bwc.R` (produces Table S1 and related plots)

---

### Data-driven simulations: GFR & AQI 🌍
- Purpose: run realistic simulations using FPCs estimated from real datasets.

- Preparation: `eigf_aqi.rds`, `eigval_aqi.rds`, `tgrid_aqi.rds` are required for AQI data simulation and `eigf_gfr.rds`, `eigval_gfr.rds` are needed for GFR data simulation. These files are generated in the real-data analysis.

- GFR (kidney transplant): `experiments/simu_gfr.R`
  - Run: `Rscript experiments/simu_gfr.R --seed 1` or use `run_reps.sh` for many seeds
  - Outputs per seed: `experiments/gfr/simu_gfr_sdXXXX.csv`, and `experiments/gfr/Theta_gfr_sdXXXX.rds` (ThetaRecord / per-pass averaged FPCs)
  - Analysis: `Rscript experiments/analyze_simu_gfr.R` → FPC plots and summary tables (`simu-gfr-fpc.pdf`)
  - Note: GFR uses `gendata.gfr()` and can depend on the private CSV (`data/gfr/GFR2023Feb01.csv`) if you want to re-generate the original estimates.

- AQI (EPA PM10): `experiments/simu_aqi.R`
  - Run: `Rscript experiments/simu_aqi.R --seed 1` or `./experiments/run_reps.sh experiments/simu_aqi.R 100`
  - Outputs per seed: `experiments/aqi/simu_aqi_sdXXXX.csv`, and `experiments/aqi/Theta_aqi_sdXXXX.rds`
  - Analysis: `Rscript experiments/analyze_simu_aqi.R` → plots such as `simu-aqi-fpc.pdf`
  - Note: AQI helper `gendata.aqi()` uses preprocessed `data/epa-aqs/aqi-us.Rda` included in the repo.
