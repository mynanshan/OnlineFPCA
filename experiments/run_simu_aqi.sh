#!/usr/bin/env bash

reps=${1:-100}

for (( i=1; i<=reps; i++ )); do
  echo "Running simulation iteration $i..."
  sbatch experiments/simu_aqi.R --seed "$i"
done
