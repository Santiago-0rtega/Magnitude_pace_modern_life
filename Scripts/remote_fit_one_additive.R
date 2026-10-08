# Fits one time-adjusted (additive) location-scale model, <moderator> + <timevar>,
# with the same specification as the primary models (remote_fit_one_primary.R):
# brms default priors with sd(es_id_model) fixed at 1, diagonal sampling
# covariance matrix, and default_mcmc_args. Supersedes run_additive.R.
# Usage (from the repository root):
#   Rscript Scripts/remote_fit_one_additive.R <moderator> <log10_years|log10_generations>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("Usage: Rscript remote_fit_one_additive.R MODERATOR TIMEVAR")
moderator <- args[[1]]
timevar   <- args[[2]]
stopifnot(moderator %in% c("disturbance", "design", "trait_type", "genphen"),
          timevar   %in% c("log10_years", "log10_generations"))

setwd(here::here())
library(here)
for (f in c("00_packages_fit", "01_paths", "05_phylogeny", "07_model_formulas",
            "08_fit_or_read_model", "10_model_diagnostics"))
  source(here::here("Scripts", paste0(f, ".R")))

dat_es   <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
A_full   <- readRDS(dir_out("phylogeny", "proceed_A_matrix.rds"))
name_map <- readRDS(dir_out("phylogeny", "proceed_name_map.rds"))
dat_es   <- apply_phylo_name_map(dat_es, name_map)

label      <- paste0(moderator, "_plus_", timevar)
model_name <- paste0(label, "_ls_additive")

dat_model <- dat_es[!is.na(dat_es[[moderator]]) & !is.na(dat_es[[timevar]]), ]
dat_model[[moderator]] <- droplevels(factor(dat_model[[moderator]]))
dat_model <- dat_model |> dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))

phylo     <- prepare_phylo_and_data(dat_model, A_full, label = label)
dat_model <- phylo$dat_model; A_mod <- phylo$A_mod; V <- phylo$V

# Main effects only, no interaction, in both the location and scale components.
formula <- build_ls_formula(paste(moderator, "+", timevar),
                            has_phylogeny = phylo$has_phylo)
priors  <- build_ls_priors(formula, dat_model, V, A = A_mod)
stopifnot(isTRUE(verify_esid_prior(priors, model_name)))

cat("=== Additive model:", label, "| n =", nrow(dat_model), "===\n")
fit <- fit_or_read_model(model_name,
  fit_fun = function() fit_ls_model(dat_model, formula, priors, V, A_mod,
                                     mcmc_args = default_mcmc_args),
  model_dir = dir_out("models"), refit = TRUE)

diag <- extract_diagnostics(fit, label, moderator)
print(diag)
if (diag$n_divergent > 0 || diag$max_rhat > 1.01 || diag$min_bulk_ess < 400)
  stop("Convergence criteria not met for ", label)
cat("CONVERGED:", label, "\n")
