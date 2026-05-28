moderator_grid <- tibble::tribble(
  ~model_id, ~moderator,           ~label,                              ~type,
  "m01",     "disturbance",        "Disturbance context",               "categorical",
  "m02",     "design",             "Comparison design",                 "categorical",
  "m03",     "log10_years",        "Elapsed time (log10 years)",        "continuous",
  "m04",     "log10_generations",  "Elapsed time (log10 generations)",  "continuous",
  "m05",     "trait_type",         "Trait type",                        "categorical",
  "m06",     "taxa",               "Taxonomic group",                   "categorical",
  "m07",     "genphen",            "Phenotypic vs genetic study",       "categorical",
  "m08",     "env_change",         "Environmental-change context",      "categorical",
  "m09",     "data_type",          "Data type",                         "categorical",
  "m10",     "transf_data",        "Transformation status",             "categorical",
  "m11",     "data_scale",         "Measurement scale",                 "categorical"
)
