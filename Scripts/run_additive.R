# Fits one time-adjusted (additive) location-scale model: <moderator> + <timevar>.
# Usage (from the repository root): Rscript Scripts/run_additive.R <moderator> <log10_years|log10_generations>
# Run from the repository root (paths resolved with here::here()).
library(here)
source(here::here("Scripts", "00_packages.R"))
source(here::here("Scripts", "01_paths.R"))
source(here::here("Scripts", "05_phylogeny.R"))
source(here::here("Scripts", "07_model_formulas.R"))
source(here::here("Scripts", "08_fit_or_read_model.R"))

args      <- commandArgs(trailingOnly = TRUE)
moderator <- args[1]   # e.g. "disturbance", "design", "trait_type", "genphen"
timevar   <- args[2]   # "log10_years" or "log10_generations"

refit_models <- FALSE

dat_es   <- readRDS(here::here("outputs", "effect_sizes", "proceed_lnm_safe.rds"))
A_full   <- readRDS(here::here("outputs", "phylogeny", "proceed_A_matrix.rds"))
name_map <- readRDS(here::here("outputs", "phylogeny", "proceed_name_map.rds"))

if (!is.null(name_map)) {
  dat_es <- apply_phylo_name_map(dat_es, name_map)
} else {
  dat_es$sp_ncbi_canonical <- dat_es$sp_ncbi
}

# --- Prepare model data ---
dat_model <- dat_es[!is.na(dat_es[[moderator]]) & !is.na(dat_es[[timevar]]), ]
dat_model <- dat_model[!is.na(dat_model$sp_ncbi_canonical), ]
dat_model[[moderator]] <- droplevels(factor(dat_model[[moderator]]))
dat_model$es_id_model  <- factor(dat_model$es_id_model)

label <- paste0(moderator, "_plus_", timevar)
cat("=== Additive model:", label, "===\n")
cat("n =", nrow(dat_model), "\n")
cat(moderator, "levels:\n"); print(table(dat_model[[moderator]]))

phylo_prep <- prepare_phylo_and_data(dat_model, A_full, label = label)
dat_model  <- phylo_prep$dat_model
A_mod      <- phylo_prep$A_mod

V <- metafor::vcalc(vi_lnM_safe, cluster = ref_id,
                    obs = es_id_model, data = dat_model, rho = 0.5)
rownames(V) <- colnames(V) <- dat_model$es_id_model

# --- Additive formula: main effects only, no interaction ---
loc_formula <- stats::as.formula(paste0(
  "yi_lnM_safe ~ ", moderator, " + ", timevar,
  " + (1 | ref_id)",
  " + (1 | gr(sp_ncbi_canonical, cov = A))",
  " + (1 | gr(es_id_model, cov = V))"
))
scl_formula <- stats::as.formula(paste0("sigma ~ ", moderator, " + ", timevar))

formula_add <- brms::bf(loc_formula, scl_formula)

priors_add <- c(
  brms::prior(normal(0, 1),   class = b),
  brms::prior(normal(0, 1),   class = b, dpar = sigma),
  brms::prior(normal(0, 0.5), class = sd),
  brms::prior(constant(1),    class = sd, group = es_id_model)
)

model_name <- paste0(label, "_ls_additive")

fit_add <- fit_or_read_model(
  model_name = model_name,
  fit_fun = function() {
    brms::brm(
      formula = formula_add,
      data    = dat_model,
      data2   = list(A = A_mod, V = V),
      prior   = priors_add,
      chains  = 4, iter = 4000, warmup = 2000, cores = 4,
      backend = "cmdstanr",
      control = list(adapt_delta = 0.97, max_treedepth = 15)
    )
  },
  refit = refit_models
)

cat("\n===", model_name, "summary ===\n")
print(summary(fit_add))
cat("\nDone.\n")
