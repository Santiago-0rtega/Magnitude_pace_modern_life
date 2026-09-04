args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript remote_fit_one_primary.R MODEL_ID")
model_id <- args[[1]]

setwd(here::here())
library(here)
for (f in c("00_packages_fit", "01_paths", "05_phylogeny", "06_model_registry",
            "07_model_formulas", "08_fit_or_read_model", "10_model_diagnostics"))
  source(here::here("chapters", "Rscripts", paste0(f, ".R")))

dat_es <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
A_full <- readRDS(dir_out("phylogeny", "proceed_A_matrix.rds"))
name_map <- readRDS(dir_out("phylogeny", "proceed_name_map.rds"))
dat_es <- apply_phylo_name_map(dat_es, name_map)

# The fittable set comes from the registry itself, so it cannot drift out of
# step with 06_model_registry.R the way a hardcoded list did. excluded_grid
# models are still fittable — m08 is rendered as an appendix chapter.
retained <- c("m00", moderator_grid$model_id, excluded_grid$model_id)
if (!model_id %in% retained)
  stop("Not a model in the registry: ", model_id,
       "\nAvailable: ", paste(retained, collapse = ", "), call. = FALSE)

if (model_id == "m00") {
  moderator <- NULL; dat_model <- dat_es; model_name <- "m00_ls_intercept_only"
} else {
  spec <- moderator_grid[moderator_grid$model_id == model_id, ]
  if (nrow(spec) == 0L) {
    # excluded_grid carries no `type` column; infer it from the data instead.
    spec <- excluded_grid[excluded_grid$model_id == model_id, ]
    spec$type <- if (is.numeric(dat_es[[spec$moderator[[1]]]]))
      "continuous" else "categorical"
  }
  stopifnot(nrow(spec) == 1L)
  moderator <- spec$moderator[[1]]
  dat_model <- dat_es[!is.na(dat_es[[moderator]]), ]
  if (spec$type[[1]] == "categorical")
    dat_model[[moderator]] <- droplevels(factor(dat_model[[moderator]]))
  model_name <- paste0(model_id, "_ls_", moderator)
}

dat_model <- dat_model |> dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
phylo <- prepare_phylo_and_data(dat_model, A_full, label = model_id)
dat_model <- phylo$dat_model; A_mod <- phylo$A_mod; V <- phylo$V
formula <- if (model_id == "m00") brms::bf(
  yi_lnM_safe ~ 1 + (1 | ref_id) + (1 | gr(sp_ncbi, cov = A)) +
    (1 | gr(es_id_model, cov = V)), sigma ~ 1
) else build_ls_formula(moderator, has_phylogeny = phylo$has_phylo)
priors <- build_ls_priors(formula, dat_model, V, A = A_mod)
stopifnot(isTRUE(verify_esid_prior(priors, model_name)))
fit <- fit_or_read_model(model_name,
  fit_fun = function() fit_ls_model(dat_model, formula, priors, V, A_mod,
                                     mcmc_args = default_mcmc_args),
  model_dir = dir_out("models"), refit = TRUE)
diag <- extract_diagnostics(fit, model_id,
  if (is.null(moderator)) "intercept" else moderator)
print(diag)
if (diag$n_divergent > 0 || diag$max_rhat > 1.01 || diag$min_bulk_ess < 400)
  stop("Convergence criteria not met for ", model_id)
cat("CONVERGED:", model_id, "\n")
