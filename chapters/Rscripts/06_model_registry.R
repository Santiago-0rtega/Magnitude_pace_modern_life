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
  "m07",     "genphen",            "Phenotypic vs genetic study",       "categorical",
  "m08",     "env_change",         "Environmental-change context",      "categorical",
  "m10",     "transf_data",        "Transformation status",             "categorical",
  "m11",     "data_scale",         "Measurement scale",                 "categorical"
)

# Moderators fit but excluded from the primary supplementary results, with
# rationale. Rendered as a documented exclusion note in ch05.
excluded_grid <- tibble::tribble(
  ~model_id, ~moderator,  ~label,             ~reason,
  "m06",     "taxa",      "Taxonomic group",
    "Every collapsing scheme tried (taxa_v2, then a 3-group Vertebrate/Invertebrate/Plant split) still fit far more slowly than every other moderator model, and never reached a usable posterior in a computationally reasonable time. The taxonomic grouping is largely redundant with the phylogenetic random effect (1 | gr(sp_ncbi, cov = A)) already in every model, so the two terms likely compete for the same signal, producing a poorly identified, slow-mixing posterior rather than a genuine model failure. Excluded as computationally infeasible given this redundancy.",
  "m09",     "data_type", "Data type",
    "The cleaned refit (m09b: index and temperature levels dropped, reference = linear) still returned 5/8000 (0.06%) divergent transitions after warmup. Although the rate is low, a divergence indicates the sampler could not fully explore the posterior in that region, so estimates from this model are not treated as reliable for reporting. Both the original (m09) and refit (m09b) fits are retained in outputs/models/ for reproducibility, but excluded from the primary results."
)
