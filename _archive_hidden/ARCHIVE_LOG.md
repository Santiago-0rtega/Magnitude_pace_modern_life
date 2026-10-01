# Archive log — `_archive_hidden/`

Generated: 2026-06-01

## Files moved to archive

| Original path | Reason |
|---|---|
| `Rscripts/Pilot_run.Rmd` | Early pilot notebook; logic superseded by `R/02–04_*.R` helpers and chapter `01_data_preparation.qmd`. Zero cross-references found in any active `.qmd` or `.R` file. |
| `Rscripts/Pilot_run.nb.html` | Rendered output of the above; no longer needed. |
| `Rscripts/SUMAMRIES.Rmd` | Alluvial-plot exploration notebook; zero cross-references found. Its alluvial summary content is not part of the final book outline. |
| `Plots/` | Empty legacy output folder; superseded by `outputs/figures/`. |
| `Rdata/` | Empty legacy model-cache folder; superseded by `outputs/models/`. |

## Ambiguous files — NOT archived

None flagged. All active R scripts (`R/00_packages.R` … `R/13_epred_draws_plots.R`) and all nine Quarto chapters are explicitly sourced or listed in `_quarto.yml`.

## Restoration

To restore any archived file, move it back to its original path relative to the repo root. No data was deleted.

## 2026-10-01 — stale artifacts (`stale_2026-10-01/`)

| Original path | Reason |
|---|---|
| `Rdata/summaries/model_summaries_emmeans.txt` | Text dump from an older fit (7,455 observations; current analysis set is 7,186). Written by `Scripts/model_summaries_emmeans.R`; not read by any `.qmd`, `.R`, or `_quarto.yml`. Moved, not deleted. |

`Rdata/tables/small_study_effects.csv` was reviewed and **kept**: it is read by
`chapters/07b_small_study_effects.qmd` as the labelled REML (metafor,
homoscedastic) cross-check, and its saved model (`Rdata/models/small_study_effects.rds`)
uses the current k = 7,186 contrasts.
