setwd(here::here())
library(here)
source(here::here("chapters", "Rscripts", "00_packages_fit.R"))
source(here::here("chapters", "Rscripts", "01_paths.R"))
source(here::here("chapters", "Rscripts", "05_phylogeny.R"))
source(here::here("chapters", "Rscripts", "06_model_registry.R"))
source(here::here("chapters", "Rscripts", "07_model_formulas.R"))
source(here::here("chapters", "Rscripts", "08_fit_or_read_model.R"))
source(here::here("chapters", "Rscripts", "09_model_summaries.R"))
source(here::here("chapters", "Rscripts", "10_model_diagnostics.R"))

refit_models <- FALSE

dat_es  <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
A_path  <- dir_out("phylogeny", "proceed_A_matrix.rds")
nm_path <- dir_out("phylogeny", "proceed_name_map.rds")

phylo_out <- list(A = readRDS(A_path), name_map = readRDS(nm_path))
A_full    <- phylo_out$A
name_map  <- phylo_out$name_map

if (!is.null(name_map)) {
  dat_es <- apply_phylo_name_map(dat_es, name_map)
} else {
  dat_es$sp_ncbi_canonical <- dat_es$sp_ncbi
}

# Collapse transf_data to untransformed vs transformed
# Drop arcsin.sqr (4 obs, too sparse)
dat_es$transf_data_v2 <- dplyr::case_when(
  dat_es$transf_data == "raw"        ~ "Untransformed",
  dat_es$transf_data == "arcsin.sqr" ~ NA_character_,   # drop (n=4)
  !is.na(dat_es$transf_data)         ~ "Transformed",
  TRUE ~ NA_character_
)

cat("transf_data_v2 counts:\n")
print(table(dat_es$transf_data_v2, useNA = "ifany"))

message("\n--- m10b: transf_data_v2 ---")

dat_model <- dat_es[!is.na(dat_es$transf_data_v2), ]
dat_model$transf_data_v2 <- droplevels(factor(dat_model$transf_data_v2))

dat_model <- dat_model |>
  dplyr::filter(!is.na(sp_ncbi_canonical)) |>
  dplyr::mutate(es_id_model = factor(es_id_model))

message("[m10b] ", nrow(dat_model), " rows after filtering")

phylo_prep <- prepare_phylo_and_data(dat_model, A_full, label = "m10b")
dat_model  <- phylo_prep$dat_model
A_mod      <- phylo_prep$A

V <- metafor::vcalc(vi_lnM_safe, cluster = ref_id,
                    obs = es_id_model, data = dat_model, rho = 0.5)
rownames(V) <- colnames(V) <- dat_model$es_id_model

formula <- build_ls_formula("transf_data_v2")
priors  <- build_ls_priors(formula, dat_model, V, A_mod)
verify_esid_prior(priors, "m10b")

fit_m10b <- fit_or_read_model(
  model_name = "m10b_ls_transf_data_v2",
  fit_fun    = function() {
    fit_ls_model(dat_model, formula, priors, V, A_mod,
                 mcmc_args = default_mcmc_args)
  },
  refit = refit_models
)

cat("\n=== m10b summary ===\n")
print(summary(fit_m10b))
message("Done.")
