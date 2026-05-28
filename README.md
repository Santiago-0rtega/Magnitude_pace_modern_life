# Magnitude pace modern life — supplementary analyses

This repository contains the Quarto book and analysis code for the lnM reanalysis of phenotypic divergence in PROCEED.

The book is the computational supplement to:

> **Human disturbance and the magnitude of phenotypic divergence** — Santiago Ortega and collaborators.

---

## Render the book

```r
quarto::quarto_render()
```

Or from the terminal:

```bash
quarto render
```

---

## Default behaviour

By default, the book **reads cached objects** from `outputs/`:

- Effect-size estimates (`outputs/effect_sizes/proceed_lnm_safe.rds`) are read without recomputing.
- Fitted brms models (`outputs/models/*.rds`) are read without refitting.

The book will render in minutes in default mode.

---

## Recompute effect sizes

Open `02_lnm_effect_sizes.qmd` and set:

```r
recompute_effect_sizes <- TRUE
```

This re-runs the SAFE lnM estimator for all contrasts. It is slow (minutes to hours depending on dataset size) and requires `orchaRd` with `orchaRd:::.safe_lnM_indep()` available.

---

## Refit models

Open `05_location_scale_model_grid.qmd` and set:

```r
refit_models <- TRUE
```

This refits all location–scale brms models from scratch using `cmdstanr`. With `chains = 4, iter = 4000`, each model takes several minutes to hours. Ensure CmdStan is installed:

```r
cmdstanr::install_cmdstan()
```

---

## Main outputs

| Path | Contents |
|---|---|
| `outputs/effect_sizes/` | SAFE lnM estimates (`.rds`) and diagnostics (`.csv`) |
| `outputs/data_clean/` | Cleaned and filtered PROCEED dataset (`.rds`) |
| `outputs/phylogeny/` | Phylogenetic tree and correlation matrix (`.rds`) |
| `outputs/models/` | Fitted brms location–scale models (`.rds`) |
| `outputs/models/sensitivity/` | Sensitivity analysis models (`.rds`) |
| `outputs/tables/` | Result and diagnostic tables (`.csv`) |
| `outputs/tables/sensitivity/` | Sensitivity analysis tables (`.csv`) |
| `outputs/figures/pdf/` | Publication-ready PDF figures |
| `outputs/figures/png/` | Web-ready PNG figures |

---

## Book structure

| Chapter | Purpose |
|---|---|
| `index.qmd` | Introduction and reproducibility guide |
| `01_data_preparation.qmd` | Data import, cleaning, inclusion criteria |
| `02_lnm_effect_sizes.qmd` | SAFE lnM computation |
| `03_descriptive_summaries.qmd` | Summary tables and distributions |
| `04_dependency_structure.qmd` | Non-independence diagnostics |
| `05_location_scale_model_grid.qmd` | Fit/read all 11 location–scale models |
| `06_location_scale_results_by_moderator.qmd` | Results for each moderator |
| `07_sensitivity_analyses.qmd` | N≥40, no-phylogeny, sys_id sensitivities |
| `08_publication_figures.qmd` | Generate and save all figures |
| `09_reproducibility.qmd` | Session info, git hash, model cache |

---

## R scripts

All reusable functions live in `R/`:

| Script | Contents |
|---|---|
| `00_packages.R` | Package loading |
| `01_paths.R` | Output directory definitions |
| `02_read_clean_data.R` | Data import and cleaning functions |
| `03_filters.R` | Inclusion criteria and derived variables |
| `04_lnm_safe.R` | SAFE lnM wrapper and caching |
| `05_phylogeny.R` | Phylogenetic scaffold via rotl |
| `06_model_registry.R` | Moderator grid (tibble) |
| `07_model_formulas.R` | brms formula and prior builders |
| `08_fit_or_read_model.R` | Cached model fitting |
| `09_model_summaries.R` | Posterior summary extraction |
| `10_model_diagnostics.R` | MCMC diagnostic extraction |
| `11_plotting_orchard_like.R` | Orchard-like plotting functions |
| `12_save_figures.R` | Dual PDF/PNG figure saving |

---

## Dependencies

```r
install.packages(c(
  "tidyverse", "here", "janitor", "brms", "cmdstanr",
  "tidybayes", "posterior", "bayesplot",
  "ape", "rotl", "patchwork", "sessioninfo"
))

# orchaRd from GitHub (provides SAFE lnM estimator)
remotes::install_github("daniel1noble/orchaRd")

# prepR4pcm — species name matching + phylogeny retrieval
# Uses a 4-stage cascade: exact → normalised → synonym (taxadb) → fuzzy
remotes::install_github("itchyshin/prepR4pcm")

# taxadb — required for the synonym-resolution step in reconcile_tree()
# Without it the pipeline falls back to exact + normalised + fuzzy matching only
install.packages("taxadb")
taxadb::td_create("col")   # download Catalogue of Life snapshot (~500 MB, once only)

# Install CmdStan
cmdstanr::install_cmdstan()
```
