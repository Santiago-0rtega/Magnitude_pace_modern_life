#!/usr/bin/env bash
set -euo pipefail
cd /home/ortegara/Documents/PACE

stamp=$(date +%Y%m%d_%H%M%S)
run_dir="logs/small_study_${stamp}"
mkdir -p "$run_dir" outputs/small_study
echo "$run_dir" > logs/current_small_study_dir.txt

models=(n_se n_v)
pids=()
for model in "${models[@]}"; do
  Rscript Scripts/remote_fit_one_small_study.R "$model" >"$run_dir/${model}.log" 2>&1 &
  pids+=("$!")
done

failed=0
for i in "${!models[@]}"; do
  if wait "${pids[$i]}"; then
    echo "OK: ${models[$i]}" | tee -a "$run_dir/status.log"
  else
    echo "FAILED: ${models[$i]}" | tee -a "$run_dir/status.log"
    failed=1
  fi
done

[[ "$failed" -eq 0 ]] || exit 1
Rscript Scripts/remote_finalize_small_study.R >"$run_dir/finalize.log" 2>&1
echo "COMPLETE: $run_dir" | tee -a "$run_dir/status.log"
