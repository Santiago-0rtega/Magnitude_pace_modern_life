#!/usr/bin/env bash
set -euo pipefail
cd /home/ortegara/Documents/PACE

stamp=$(date +%Y%m%d_%H%M%S)
run_dir="logs/small_study_brms_${stamp}"
mkdir -p "$run_dir" outputs/tables/small_study_brms
echo "$run_dir" > logs/current_small_study_brms_dir.txt

models=(n_se_sigma1 n_se_sigmax n_v_sigma1 n_v_sigmax)
pids=()
for model in "${models[@]}"; do
  Rscript Scripts/remote_fit_one_small_study_brms.R "$model" >"$run_dir/${model}.log" 2>&1 &
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

if [[ "$failed" -ne 0 ]]; then
  echo "FAILED: one or more brms models" | tee -a "$run_dir/status.log"
  exit 1
fi

echo "COMPLETE: $run_dir" | tee -a "$run_dir/status.log"
