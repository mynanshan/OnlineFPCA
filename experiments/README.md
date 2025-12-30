# Simulation Studies

Experiment files:

* `fpca1d.R`: main 1D simulation
* `fpca2d.R`: main 2D simulation
* `ci-simu1d.R`: 1D simulation with pointwise CI
* `ci-simu2d.R`: 2D simulation with pointwise CI
* `check_abv.R`, `plot_abv.R`: Check how $\omega$ and block size affects 1D simulation performance
* `check_bwc.R`: Check how $(C,W,B)$ affects 1D simulation performance
* `check_nbatch.R`: Check how batch size affects 1D simulation performance
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