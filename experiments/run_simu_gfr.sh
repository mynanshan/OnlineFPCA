#!/usr/bin/env bash

reps=${1:-100}

for (( i=1; i<=reps; i++ )); do
  echo "Running simulation iteration $i..."
  sbatch examples/simu_simple.R \
    --seed        "$i"
done