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

# ── Scheme C: collapse taxa into 3 broad, well-populated groups ───────────────
# Vertebrate = Fish+Bird+Mammal+Reptile+Amphibian, Invertebrate = Arthropod+
# Mollusc+Annelid, Plant. Reference level = Vertebrate (largest group).
dat_es$taxa3 <- dplyr::case_when(
  dat_es$taxa %in% c("Fish", "Bird", "Mammal", "Reptile", "Amphibian") ~ "Vertebrate",
  dat_es$taxa %in% c("Arthropod", "Mollusc", "Annelid")                ~ "Invertebrate",
  dat_es$taxa == "Plant"                                               ~ "Plant",
  TRUE ~ NA_character_
)
dat_es$taxa3 <- factor(dat_es$taxa3,
                       levels = c("Vertebrate", "Invertebrate", "Plant"))

cat("taxa3 counts:\n"); print(table(dat_es$taxa3, useNA = "ifany"))

# ── Generic fitter (identical to run_m06b.R) ─────────────────────────────────
fit_new_model <- function(mod_id, moderator, mod_label) {
  message("\n--- ", mod_id, ": ", moderator, " ---")

  dat_model <- dat_es[!is.na(dat_es[[moderator]]), ]
  dat_model[[moderator]] <- droplevels(factor(dat_model[[moderator]]))

  sparse <- names(which(table(dat_model[[moderator]]) < 5))
  if (length(sparse)) {
    message("Sparse levels (<5 obs): ", paste(sparse, collapse = ", "))
    dat_model <- dat_model[!dat_model[[moderator]] %in% sparse, ]
    dat_model[[moderator]] <- droplevels(dat_model[[moderator]])
  }

  dat_model <- dat_model |>
    dplyr::filter(!is.na(sp_ncbi_canonical)) |>
    dplyr::mutate(es_id_model = factor(es_id_model))

  message("[", mod_id, "] ", nrow(dat_model), " rows after filtering")

  phylo_prep <- prepare_phylo_and_data(dat_model, A_full, label = mod_id)
  dat_model  <- phylo_prep$dat_model
  A_mod      <- phylo_prep$A

  V <- metafor::vcalc(vi_lnM_safe, cluster = ref_id,
                      obs = es_id_model, data = dat_model, rho = 0.5)
  rownames(V) <- colnames(V) <- dat_model$es_id_model

  formula <- build_ls_formula(moderator)
  priors  <- build_ls_priors(formula, dat_model, V, A_mod)
  verify_esid_prior(priors, mod_id)

  fit_or_read_model(
    model_name = paste0(mod_id, "_ls_", moderator),
    fit_fun    = function() {
      fit_ls_model(dat_model, formula, priors, V, A_mod,
                   mcmc_args = default_mcmc_args)
    },
    refit = refit_models
  )
}

# ── Fit m06c (taxa3) ─────────────────────────────────────────────────────────
fit_m06c <- fit_new_model("m06c", "taxa3", "Taxonomic group (broad)")
cat("\n=== m06c summary ===\n"); print(summary(fit_m06c))

message("\nDone.")
