# Primary moderator grid: models included in the reported supplementary results.
# m06 (taxa), m08 (env_change), m09 (data_type), m10 (transf_data) and m11
# (data_scale) were fit but are EXCLUDED here — see
# excluded_grid below for rationale. Keeping them out of this table means they
# are never loaded into the grid status / effect-summary tables that drive the
# publication figures (ch08).
moderator_grid <- tibble::tribble(
  ~model_id, ~moderator,           ~label,                              ~type,
  "m01",     "disturbance",        "Disturbance context",               "categorical",
  "m02",     "design",             "Comparison design",                 "categorical",
  "m03",     "log10_years",        "Elapsed time (log10 years)",        "continuous",
  "m04",     "log10_generations",  "Elapsed time (log10 generations)",  "continuous",
  "m05",     "trait_type",         "Trait type",                        "categorical",
  "m07",     "genphen",            "Phenotypic vs genetic study",       "categorical"
)

# Non-primary model families. This table is rendered once, at the end of ch05.
excluded_grid <- tibble::tribble(
  ~model_id, ~moderator,  ~label,             ~reason,
  "m06",     "taxa",      "Taxonomic group",
    "The original fit completed only 2 of 4 chains (max Rhat 1.05; minimum sigma ESS 59). Two collapsed-category refits were attempted, but neither produced a usable converged posterior within the allocated runtime.",
  "m08",     "env_change", "Environmental-change context",
    "The corrected-data fit had max Rhat 1.0143, minimum bulk ESS 780, minimum tail ESS 663, and 0 divergences. It did not meet the pre-specified Rhat limit.",
  "m09",     "data_type", "Data type",
    "The original fit had bulk ESS near 350 for sigma terms. The simplified m09b refit removed two sparse levels but retained 5 of 8000 post-warmup divergent transitions.",
  "m10",     "transf_data", "Transformation status",
    "The original seven-level fit stopped because of memory limits. The simplified two-level m10b refit had max Rhat 1.0128, minimum bulk ESS 230, minimum tail ESS 363, and 0 divergences.",
  "m11",     "data_scale", "Measurement scale",
    "Dropped by the authors; not part of the reported model set. The fitted model file is retained but its results are not reported."
)
