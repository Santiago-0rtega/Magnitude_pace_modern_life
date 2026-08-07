args <- commandArgs(trailingOnly = TRUE)
valid_models <- c("n_se_sigma1", "n_se_sigmax", "n_v_sigma1", "n_v_sigmax")
if (length(args) != 1L || !args[[1]] %in% valid_models) {
  stop("Usage: Rscript remote_fit_one_small_study_brms.R ",
       paste(valid_models, collapse = "|"))
}
model_id <- args[[1]]
extended_refit <- identical(tolower(Sys.getenv("EXTENDED_REFIT", "false")), "true")
artifact_suffix <- Sys.getenv("ARTIFACT_SUFFIX", "")
artifact_id <- paste0(model_id, artifact_suffix)

setwd("/home/ortegara/Documents/PACE")
library(here)
for (f in c("00_packages", "01_paths", "05_phylogeny", "07_model_formulas",
            "08_fit_or_read_model", "10_model_diagnostics")) {
  source(here::here("R", paste0(f, ".R")))
}

predictor <- if (startsWith(model_id, "n_se")) "n_se" else "n_v"
sigma_predictor <- endsWith(model_id, "sigmax")
model_name <- paste0("small_study_brms_", artifact_id)

dat_es <- readRDS(here::here("outputs", "effect_sizes", "proceed_lnm_safe.rds"))
A_full <- readRDS(here::here("outputs", "phylogeny", "proceed_A_matrix.rds"))
name_map <- readRDS(here::here("outputs", "phylogeny", "proceed_name_map.rds"))
dat_es <- apply_phylo_name_map(dat_es, name_map) |>
  dplyr::filter(
    is.finite(yi_lnM_safe), is.finite(vi_lnM_safe), vi_lnM_safe > 0,
    is.finite(n1), n1 > 0, is.finite(n2), n2 > 0
  ) |>
  dplyr::mutate(
    n0 = (n1 * n2) / (n1 + n2),
    n_se = 1 / sqrt(n0),
    n_v = 1 / n0,
    es_id_model = factor(seq_len(dplyr::n()))
  )

phylo <- prepare_phylo_and_data(dat_es, A_full, label = model_id)
dat_model <- phylo$dat_model
A_mod <- phylo$A_mod
V <- phylo$V
stopifnot(phylo$has_phylo, nrow(dat_model) > 0L)

location_formula <- stats::as.formula(paste0(
  "yi_lnM_safe ~ ", predictor,
  " + (1 | ref_id)",
  " + (1 | gr(sp_ncbi, cov = A))",
  " + (1 | gr(es_id_model, cov = V))"
))
sigma_formula <- if (sigma_predictor) {
  stats::as.formula(paste("sigma ~", predictor))
} else {
  sigma ~ 1
}
formula <- brms::bf(location_formula, sigma_formula)

cat("MODEL:", model_id, "\n")
cat("ARTIFACT:", artifact_id, "\n")
cat("N:", nrow(dat_model), "\n")
cat("LOCATION:", paste(deparse(location_formula), collapse = " "), "\n")
cat("SCALE:", paste(deparse(sigma_formula), collapse = " "), "\n")

priors <- build_ls_priors(formula, dat_model, V, A = A_mod)
stopifnot(isTRUE(verify_esid_prior(priors, model_name)))
print(priors)

mcmc_args <- if (extended_refit) {
  list(
    chains = 4,
    cores = 4,
    iter = 8000,
    warmup = 4000,
    seed = 124,
    backend = "cmdstanr",
    control = list(adapt_delta = 0.995, max_treedepth = 18)
  )
} else {
  default_mcmc_args
}
cat("MCMC SETTINGS:\n")
print(mcmc_args)

fit <- fit_or_read_model(
  model_name,
  fit_fun = function() {
    fit_ls_model(
      dat_model, formula, priors, V, A_mod,
      mcmc_args = mcmc_args
    )
  },
  model_dir = here::here("outputs", "models"),
  refit = TRUE
)

diag <- extract_diagnostics(
  fit, artifact_id, predictor,
  max_treedepth = mcmc_args$control$max_treedepth
)
out_dir <- here::here("outputs", "tables", "small_study_brms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(diag, file.path(out_dir, paste0(artifact_id, "_diagnostics.csv")))

fe <- as.data.frame(brms::fixef(fit))
fe$term <- rownames(fe)
fe$model_id <- artifact_id
fe <- fe[, c("model_id", "term", setdiff(names(fe), c("model_id", "term")))]
readr::write_csv(fe, file.path(out_dir, paste0(artifact_id, "_fixed_effects.csv")))

print(diag)
print(summary(fit))
if (diag$n_divergent > 0 || diag$max_rhat > 1.01 ||
    diag$min_bulk_ess < 400 || diag$min_tail_ess < 400) {
  stop("Convergence criteria not met for ", artifact_id)
}
cat("CONVERGED:", artifact_id, "\n")
