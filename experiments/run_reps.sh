#!/usr/bin/env bash
set -euo pipefail

module load r/4.5.0

# Usage
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 path/to/script.R [n_reps] [-- script-args...]"
  exit 1
fi

script=$1; shift

# sanity check: does the R script exist?
if [[ ! -f "$script" ]]; then
  echo "Error: file '$script' not found."
  exit 1
fi

# default reps
reps=100

# If next arg is an integer, treat it as n_reps and shift it off
if [[ $# -gt 0 && "$1" =~ ^[0-9]+$ ]]; then
  reps=$1
  shift
fi

# Everything remaining will be forwarded to the R script
forward_args=("$@")

for (( seed=1; seed<=reps; seed++ )); do
  echo "Submitting simulation iteration $seed..."
  ## If using a cluster with `sbatch`, use the line below to run remotely:
  # sbatch -- "$(realpath "$script")" --seed "$seed" "${forward_args[@]}"
  ## If not on a cluster with `sbatch`, use the line below to run locally:
  Rscript "$(realpath "$script")" --seed "$seed" "${forward_args[@]}"
done
