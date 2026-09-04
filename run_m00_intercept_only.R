setwd(here::here())
library(here)
source(here::here("chapters", "Rscripts", "00_packages_fit.R"))
source(here::here("chapters", "Rscripts", "01_paths.R"))
source(here::here("chapters", "Rscripts", "05_phylogeny.R"))
source(here::here("chapters", "Rscripts", "07_model_formulas.R"))
source(here::here("chapters", "Rscripts", "08_fit_or_read_model.R"))

dat_es    <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
A_path    <- dir_out("phylogeny", "proceed_A_matrix.rds")
nm_path   <- dir_out("phylogeny", "proceed_name_map.rds")

phylo_out <- list(A = readRDS(A_path), name_map = readRDS(nm_path))
A_full    <- phylo_out$A
name_map  <- phylo_out$name_map

if (!is.null(name_map)) {
  dat_es <- apply_phylo_name_map(dat_es, name_map)
} else {
  dat_es$sp_ncbi_canonical <- dat_es$sp_ncbi
}

# No moderator filter — use all rows
dat_model <- dat_es
cat("Total rows:", nrow(dat_model), "\n")

# Run through the same phylo + V/A preparation as other models
phylo_prep <- prepare_phylo_and_data(dat_model, A_full, label = "m00")
dat_model  <- phylo_prep$dat_model
A_mod      <- phylo_prep$A
cat("Rows after phylo prep:", nrow(dat_model), "\n")

V <- metafor::vcalc(vi_lnM_safe, cluster = ref_id,
                    obs = es_id_model, data = dat_model, rho = 0.5)
rownames(V) <- colnames(V) <- dat_model$es_id_model

formula_m00 <- brms::bf(
  yi_lnM_safe ~ 1 + (1 | ref_id) +
    (1 | gr(sp_ncbi, cov = A)) +
    (1 | gr(es_id_model, cov = V)),
  sigma ~ 1
)

priors <- build_ls_priors(formula_m00, dat_model, V, A_mod)
verify_esid_prior(priors, "m00")

fit_m00 <- fit_or_read_model(
  model_name = "m00_ls_intercept_only",
  fit_fun    = function() {
    fit_ls_model(dat_model, formula_m00, priors, V, A_mod,
                 mcmc_args = default_mcmc_args)
  }
)

cat("\n=== m00 intercept-only summary ===\n")
print(summary(fit_m00))
message("Done.")
