#!/usr/bin/env bash
set -euo pipefail
cd /home/anonymous/Documents/PACE
stamp=$(date +%Y%m%d_%H%M%S)
run_dir="logs/rebuild_${stamp}"
mkdir -p "$run_dir" outputs/tables/sensitivity/parts
echo "$run_dir" > logs/current_rebuild_dir.txt

Rscript Scripts/remote_rebuild_safe.R >"$run_dir/rebuild_safe.log" 2>&1

models=(m00 m01 m02 m03 m04 m05 m07 m08 m11)
pids=()
for model in "${models[@]}"; do
  Rscript Scripts/remote_fit_one_primary.R "$model" >"$run_dir/${model}.log" 2>&1 &
  pids+=("$!")
done
failed=0
for i in "${!models[@]}"; do
  if wait "${pids[$i]}"; then echo "OK primary: ${models[$i]}" | tee -a "$run_dir/status.log"
  else echo "FAILED primary: ${models[$i]}" | tee -a "$run_dir/status.log"; failed=1; fi
done
[[ "$failed" -eq 0 ]] || exit 1

# Refit all 21 established sensitivity models against the corrected dataset.
rm -f outputs/tables/sensitivity/parts/*.csv
CAP=9 REFIT_SENSITIVITY=true bash run_sensitivity_parallel.sh >"$run_dir/sensitivity.log" 2>&1

Rscript Scripts/precompute_grid_summaries.R >"$run_dir/precompute_grid.log" 2>&1
Rscript Scripts/precompute_heterogeneity.R >"$run_dir/precompute_heterogeneity.log" 2>&1
Rscript Scripts/precompute_emmeans_contrasts.R >"$run_dir/precompute_contrasts.log" 2>&1
FORCE=1 Rscript Scripts/precompute_epred_draws.R >"$run_dir/precompute_epred.log" 2>&1
ORCHARD_MODEL_IDS=m00,m01,m02,m03,m04,m05,m07,m08,m11 \
  Rscript Scripts/plot_orchard_style_v2.R >"$run_dir/plot_orchard.log" 2>&1
echo "COMPLETE: $run_dir" | tee -a "$run_dir/status.log"
