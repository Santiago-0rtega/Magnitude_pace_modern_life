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

Open `chapters/02_lnm_effect_sizes.qmd` and set:

```r
recompute_effect_sizes <- TRUE
```

This re-runs the SAFE lnM estimator for all contrasts. It is slow (minutes to hours depending on dataset size) and requires `orchaRd` with `orchaRd:::.safe_lnM_indep()` available.

---

## Refit models

Open `chapters/05_location_scale_model_grid.qmd` and set:

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
| `chapters/` | Quarto source chapters, including `chapters/models/` for model-specific chapters |
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
| `chapters/index.qmd` | Introduction and reproducibility guide |
| `chapters/01_data_preparation.qmd` | Data import, cleaning, inclusion criteria |
| `chapters/02_lnm_effect_sizes.qmd` | SAFE lnM computation |
| `chapters/03_descriptive_summaries.qmd` | Summary tables and distributions |
| `chapters/04_dependency_structure.qmd` | Non-independence diagnostics |
| `chapters/04b_phylogeny.qmd` | Phylogenetic tree construction and species-name reconciliation |
| `chapters/05_location_scale_model_grid.qmd` | Fit/read location-scale models |
| `chapters/models/m00.qmd` to `chapters/models/m13.qmd` | Model-specific results chapters; dropped `m06` and `m09` are kept as source files but omitted from `_quarto.yml` |
| `chapters/07_sensitivity_analyses.qmd` | N>=40, no-phylogeny, sys_id sensitivities |
| `chapters/08_publication_figures.qmd` | Generate and save all figures |
| `chapters/09_reproducibility.qmd` | Session info, git hash, model cache |

---

## Interpretation of lnM

### What lnM measures

lnM is the **log magnitude of phenotypic separation** between two groups. It measures how far apart two group means are relative to their typical within-group variation, without retaining the sign of the trait difference. It is therefore a magnitude effect size, not a directional contrast between a focal and reference population.

For a two-group comparison:

$$
\ln M = \ln\left(\frac{s_B}{s_W}\right),
$$

where `s_B` is the between-group separation component and `s_W` is the pooled within-group standard deviation.

Because lnM is on a log ratio scale, it can take negative values (separation is less than one pooled within-group SD) and positive values (separation is greater than one pooled within-group SD). The natural reference point is lnM = 0, i.e. between-group separation of about one pooled within-group SD.

The sign of lnM is **not** the direction of trait change. Human disturbance can increase or decrease a trait mean depending on system and trait; lnM asks about the magnitude of separation, not whether the first group is larger than the second.

### Location submodel

The location submodel estimates the **mean lnM** for each moderator level (posterior mean of the fixed effect). Interpretation:

- **Intercept**: mean lnM for the reference category of the moderator, after accounting for study (ref_id), phylogeny (sp_ncbi), and sampling variance (es_id_model).
- **Coefficient for level k**: difference in mean lnM between level k and the reference category. Positive = greater magnitude of separation in level k; negative = smaller magnitude of separation.
- **Overall baseline (m00 intercept ≈ −0.06, 95% CrI [−0.97, 0.84])**: on average across all contexts, the magnitude of phenotypic divergence is indistinguishable from zero on the log scale (i.e., M ≈ 1 SD). There is no consistent directional signal of accelerated or decelerated divergence.

### Scale submodel

The scale submodel estimates **log residual SD (log σ)** — the within-group heterogeneity in lnM after accounting for the location fixed effects and random effects. Interpretation:

- **sigma_Intercept**: log residual SD for the reference category.
- **Coefficient for level k**: difference in log SD between level k and the reference. Positive = more heterogeneous responses (studies within that level disagree more); negative = more consistent responses.
- **Overall baseline (m00 sigma_Intercept ≈ −1.09, 95% CrI [−1.11, −1.06])**: residual SD ≈ exp(−1.09) ≈ 0.34 lnM units — substantial unexplained variance remains after accounting for study and phylogeny.

The scale submodel is often more informative than the location: even when the mean lnM does not differ between groups, the *consistency* of divergence may differ substantially.

---

## Model decisions

| Model | Moderator | Type | Decision / notes |
|---|---|---|---|
| **m00** | — | Intercept-only | Baseline model. Overall mean lnM ≈ 0 (inconclusive); sigma_Intercept ≈ −1.09 (credibly negative). Study and phylogeny random effects both ~0.8 SD. |
| **m01** | `disturbance` | Categorical | All 7 original categories retained. Reference = Climate change. |
| **m02** | `design` | Categorical | Allochronic (same population, time series) vs Synchronic (diverged populations). Reference = Allochronic. |
| **m03** | `log10_years` | Continuous | Elapsed time in log₁₀ years. Strong positive location slope (+0.31) and negative scale slope (−0.16): longer studies find larger but more consistent divergence. |
| **m04** | `log10_generations` | Continuous | Elapsed time in log₁₀ generations. Same pattern as m03 (slope +0.29, scale −0.15). |
| **m05** | `trait_type` | Categorical | All 8 original types retained. Reference = Behaviour. |
| **m06** | `taxa` | Categorical | **Replaced by m06b.** Original 9-level model failed: only 2/4 chains completed, Rhat up to 1.05, ESS as low as 59 for sigma parameters. Root cause: Amphibian (reference, n=23) had near-zero within-group variance, causing sigma to collapse to −∞ and creating an unidentifiable funnel geometry. |
| **m06b** | `taxa_v2` | Categorical | Collapsed to 6 groups using `v5_taxa_fine`: **Fish** (Actinopterygii, n=3,120), **Plant** (Streptophyta, n=1,772), **Insect** (Insecta, n=1,036), **Bird** (Aves, n=847), **Mammal** (Mammalia, n=567), **Reptile** (Reptilia, n=92). Dropped Crustacea (29), Amphibia (23), Gastropoda (13), Mollusca (12) — 77 obs (~1% of data) — as insufficiently powered for the sigma submodel. Reference = Bird. |
| **m07** | `genphen` | Categorical | Genetic (common garden / QG) vs Phenotypic (wild-measured). Reference = Genetic. |
| **m08** | `env_change` | Categorical | Novel (defined start point) vs Ongoing (measured within). Reference = Novel. |
| **m09** | `data_type` | Categorical | **Replaced by m09b.** Original 12-level model had low ESS for all sigma parameters (~350 Bulk) and a problematic reference (ad_ratio, n=83, wide sigma CI). index (n=14) and temperature (n=19) too sparse for sigma submodel (same geometry as m06 failure). |
| **m09b** | `data_type_v2` | Categorical | Dropped index (n=14) and temperature (n=19). Reference changed to **linear** (n=4,148). 10 levels retained: linear, cube (3D), count, proportion, rate, time, other, date, ad_ratio, area (2D). |
| **m10** | `transf_data` | Categorical | **Collapsed to 2 levels (m10b); sensitivity analysis only.** Original 7-level model failed (OOM). Collapsed to **Untransformed** (raw, n=7,134) vs **Transformed** (n=354). Extreme imbalance (~95% raw) leaves the sigma posterior for Transformed poorly anchored (Bulk ESS ~230). Not included as a main result. |
| **m11** | `data_scale` | Categorical | Ratio scale (true zero) vs Interval scale (arbitrary zero). Reference = Interval. |
| **m12** | `log10_years × disturbance` | Interaction | Calendar-time slope allowed to vary by disturbance type (location + scale). All 7 disturbance levels retained. Note: Climate change (n=240, log₁₀ year range=1.29) and Response to introductions (n=296, range=1.18) have narrow year coverage — slopes for these levels will be less precise. If convergence fails, simplified version (collapsed disturbance) moved to sensitivity analysis. Reference = Climate change. |
| **m13** | `log10_generations × disturbance` | Interaction | Evolutionary-time slope allowed to vary by disturbance type (location + scale). Same structure as m12 but in generational time. All 7 disturbance levels retained with same caveats. Reference = Climate change. |

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
