# remote_fit_leave_oke_out.R
#
# Protocol-deviation sensitivity fit: refits m01 (disturbance-context
# location-scale model) after excluding every Oke et al. 2020 contrast
# (ref_id == "p209"), using the SAME formula, priors, and MCMC settings as
# the primary m01 fit (chains/iter/warmup/adapt_delta/max_treedepth all match
# default_mcmc_args). Run on remote-server:
#
#   Rscript remote_fit_leave_oke_out.R
#
# Saves outputs/models/sensitivity/m01_ls_disturbance_leave_oke_out.rds

setwd("/home/anonymous/Documents/PACE")
library(here)
for (f in c("00_packages", "01_paths", "05_phylogeny", "06_model_registry",
            "07_model_formulas", "08_fit_or_read_model", "09_model_summaries",
            "10_model_diagnostics"))
  source(here::here("Scripts", paste0(f, ".R")))

dat_es   <- readRDS(here::here("outputs", "effect_sizes", "proceed_lnm_safe.rds"))
A_full   <- readRDS(here::here("outputs", "phylogeny", "proceed_A_matrix.rds"))
name_map <- readRDS(here::here("outputs", "phylogeny", "proceed_name_map.rds"))
dat_es <- apply_phylo_name_map(dat_es, name_map)

moderator  <- "disturbance"
model_name <- "m01_ls_disturbance_leave_oke_out"

n_before_excl <- sum(!is.na(dat_es[[moderator]]))
dat_model <- dat_es[!is.na(dat_es[[moderator]]), ]
dat_model <- dat_model |> dplyr::filter(ref_id != "p209")   # exclude Oke et al. 2020
n_after_excl <- nrow(dat_model)
message("Rows before Oke exclusion: ", n_before_excl)
message("Rows after Oke exclusion : ", n_after_excl,
        "  (removed ", n_before_excl - n_after_excl, ")")

dat_model[[moderator]] <- droplevels(factor(dat_model[[moderator]]))
message("Disturbance level counts (no-Oke dataset):")
print(table(dat_model[[moderator]]))

dat_model <- dat_model |> dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
phylo <- prepare_phylo_and_data(dat_model, A_full, label = "m01_leave_oke_out")
dat_model <- phylo$dat_model; A_mod <- phylo$A_mod; V <- phylo$V

message("Final n rows used in fit: ", nrow(dat_model))
message("Final n species: ", dplyr::n_distinct(dat_model$sp_ncbi))
message("has_phylo: ", phylo$has_phylo)

formula <- build_ls_formula(moderator, has_phylogeny = phylo$has_phylo)
priors  <- build_ls_priors(formula, dat_model, V, A = A_mod)
stopifnot(isTRUE(verify_esid_prior(priors, model_name)))

t0 <- Sys.time()
fit <- fit_or_read_model(model_name,
  fit_fun = function() fit_ls_model(dat_model, formula, priors, V, A_mod,
                                     mcmc_args = default_mcmc_args),
  model_dir = here::here("outputs", "models", "sensitivity"), refit = TRUE)
message("Fit wall time: ", format(Sys.time() - t0))

dir.create(here::here("outputs", "tables", "sensitivity"),
          recursive = TRUE, showWarnings = FALSE)

diag <- extract_diagnostics(fit, "m01_leave_oke_out", moderator)
print(diag)
readr::write_csv(diag, here::here("outputs", "tables", "sensitivity",
                                  "m01_leave_oke_out_diag.csv"))

fe <- extract_fixed_effects(fit, "m01_leave_oke_out", moderator)
readr::write_csv(fe, here::here("outputs", "tables", "sensitivity",
                                "m01_leave_oke_out_fe.csv"))

if (diag$n_divergent > 0 || diag$max_rhat > 1.01 || diag$min_bulk_ess < 400) {
  message("NOTE: convergence criteria not met for m01_leave_oke_out (see diagnostics above).")
} else {
  message("CONVERGED: m01_leave_oke_out")
}
cat("DONE: m01_leave_oke_out\n")
