#!/usr/bin/env bash

# usage check
if [[ -z "$1" ]]; then
  echo "Usage: $0 path/to/script.R [n_reps]"
  exit 1
fi

script="$1"
reps="${2:-100}"

# sanity check: does the R script exist?
if [[ ! -f "$script" ]]; then
  echo "Error: file '$script' not found."
  exit 1
fi

for (( seed=1; seed<=reps; seed++ )); do
  echo "Submitting simulation iteration $seed..."
  sbatch $(realpath "$script") --seed $seed
done
