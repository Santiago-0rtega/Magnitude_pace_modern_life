# Primary moderator grid: models included in the reported supplementary results.
# m06 (taxa) and m09 (data_type) were fit but are EXCLUDED here — see
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

# Moderators fit but excluded from the primary supplementary results, with
# rationale. Rendered as a documented exclusion note in ch05.
excluded_grid <- tibble::tribble(
  ~model_id, ~moderator,  ~label,             ~reason,
  "m06",     "taxa",      "Taxonomic group",
    "The original taxonomic fit and two collapsed-category refits did not complete within the allocated runtime. No converged posterior was available, so this moderator is excluded from reported estimates.",
  "m08",     "env_change", "Environmental-change context",
    "The corrected-data fit had max Rhat = 1.0143 and did not meet the pre-specified convergence criterion. It is excluded from primary results and retained only as a diagnostic appendix.",
  "m09",     "data_type", "Data type",
    "The cleaned refit returned 5/8000 divergent transitions after warmup and is excluded from reported estimates.",
  "m10",     "transf_data", "Transformation status",
    "The cleaned refit had max Rhat = 1.0128 and minimum bulk ESS = 230, so it did not meet the convergence criteria and is excluded from the reported results.",
  "m11",     "data_scale", "Measurement scale",
    "Measurement scale is a property of the source measurement units rather than a biological moderator. Its converged fit is retained as a sensitivity analysis of whether lnM differs between interval- and ratio-scale traits."
)
