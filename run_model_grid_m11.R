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

# Start from m07 only
moderator_grid <- moderator_grid[moderator_grid$model_id == "m11", ]

dat_es <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
cat("Total contrasts:", nrow(dat_es), "\n")

A_path        <- dir_out("phylogeny", "proceed_A_matrix.rds")
name_map_path <- dir_out("phylogeny", "proceed_name_map.rds")
phylo_out <- list(A = readRDS(A_path), name_map = readRDS(name_map_path))
A_full    <- phylo_out$A
name_map  <- phylo_out$name_map

if (!is.null(name_map)) {
  dat_es <- apply_phylo_name_map(dat_es, name_map)
} else {
  dat_es$sp_ncbi_canonical <- dat_es$sp_ncbi
}

all_fe      <- list()
all_diag    <- list()
status_rows <- list()

for (i in seq_len(nrow(moderator_grid))) {
  row        <- moderator_grid[i, ]
  mod_id     <- row$model_id
  moderator  <- row$moderator
  mod_label  <- row$label
  mod_type   <- row$type
  model_file <- paste0(mod_id, "_ls_", moderator)

  message("\n--- ", mod_id, ": ", moderator, " ---")

  if (!moderator %in% names(dat_es)) {
    status_rows[[i]] <- tibble::tibble(model_id = mod_id, moderator = moderator,
      label = mod_label, type = mod_type, n_rows = 0L, n_levels = NA_integer_,
      n_missing_excluded = NA_integer_, model_path = NA_character_,
      status = "skipped: variable absent", notes = "Moderator not found in dataset.",
      max_rhat = NA_real_, min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
      n_divergent = NA_integer_, max_treedepth_hits = NA_integer_)
    next
  }

  n_before     <- nrow(dat_es)
  dat_model    <- dat_es[!is.na(dat_es[[moderator]]), ]
  n_missing_ex <- n_before - nrow(dat_model)

  if (mod_type == "categorical") {
    dat_model[[moderator]] <- droplevels(factor(dat_model[[moderator]]))
  }

  n_levels <- NA_integer_
  if (mod_type == "categorical") {
    level_counts <- table(dat_model[[moderator]])
    n_levels     <- length(level_counts)
    sparse       <- level_counts[level_counts < 5]
    if (length(sparse) > 0)
      message("Sparse levels (<5 obs) in ", moderator, ": ", paste(names(sparse), collapse = ", "))
    if (all(level_counts < 3)) {
      status_rows[[i]] <- tibble::tibble(model_id = mod_id, moderator = moderator,
        label = mod_label, type = mod_type, n_rows = nrow(dat_model), n_levels = n_levels,
        n_missing_excluded = n_missing_ex, model_path = NA_character_,
        status = "skipped: too sparse", notes = "All levels have fewer than 3 observations.",
        max_rhat = NA_real_, min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
        n_divergent = NA_integer_, max_treedepth_hits = NA_integer_)
      next
    }
  }

  if (nrow(dat_model) < 10) {
    status_rows[[i]] <- tibble::tibble(model_id = mod_id, moderator = moderator,
      label = mod_label, type = mod_type, n_rows = nrow(dat_model), n_levels = n_levels,
      n_missing_excluded = n_missing_ex, model_path = NA_character_,
      status = "skipped: insufficient data", notes = paste("Only", nrow(dat_model), "rows."),
      max_rhat = NA_real_, min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
      n_divergent = NA_integer_, max_treedepth_hits = NA_integer_)
    next
  }

  dat_model <- dat_model |>
    dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))

  phylo     <- prepare_phylo_and_data(dat_model, A_full, label = mod_id)
  dat_model <- phylo$dat_model
  A_mod     <- phylo$A_mod
  V         <- phylo$V
  has_phylo <- phylo$has_phylo

  formula <- build_ls_formula(moderator, has_phylogeny = has_phylo)
  priors  <- build_ls_priors(formula, dat_model, V, A = A_mod)

  prior_ok <- verify_esid_prior(priors, model_name = mod_id)
  if (!prior_ok) {
    status_rows[[i]] <- tibble::tibble(model_id = mod_id, moderator = moderator,
      label = mod_label, type = mod_type, n_rows = nrow(dat_model), n_levels = n_levels,
      n_missing_excluded = n_missing_ex, model_path = NA_character_,
      status = "skipped: prior check failed", notes = "constant(1) not confirmed.",
      max_rhat = NA_real_, min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
      n_divergent = NA_integer_, max_treedepth_hits = NA_integer_)
    next
  }

  model_path_full <- file.path(dir_out("models"), paste0(model_file, ".rds"))

  fit <- tryCatch({
    fit_or_read_model(
      model_name = model_file,
      fit_fun    = function() {
        fit_ls_model(dat_model, formula, priors, V = V, A = A_mod,
                     mcmc_args = default_mcmc_args)
      },
      model_dir = dir_out("models"),
      refit     = refit_models
    )
  }, error = function(e) {
    message("Model fitting failed for ", moderator, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(fit)) {
    status_rows[[i]] <- tibble::tibble(model_id = mod_id, moderator = moderator,
      label = mod_label, type = mod_type, n_rows = nrow(dat_model), n_levels = n_levels,
      n_missing_excluded = n_missing_ex, model_path = model_path_full,
      status = "failed", notes = "Model returned NULL.",
      max_rhat = NA_real_, min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
      n_divergent = NA_integer_, max_treedepth_hits = NA_integer_)
    next
  }

  fe   <- extract_fixed_effects(fit, mod_id, moderator)
  diag <- tryCatch(
    extract_diagnostics(fit, mod_id, moderator),
    error = function(e) {
      tibble::tibble(model_id = mod_id, moderator = moderator,
        max_rhat = NA_real_, min_bulk_ess = NA_real_,
        min_tail_ess = NA_real_, n_divergent = NA_integer_,
        max_treedepth_hits = NA_integer_)
    }
  )

  all_fe[[mod_id]]   <- split_location_scale(fe)
  all_diag[[mod_id]] <- diag

  status_rows[[i]] <- tibble::tibble(
    model_id = mod_id, moderator = moderator, label = mod_label,
    type = mod_type, n_rows = nrow(dat_model), n_levels = n_levels,
    n_missing_excluded = n_missing_ex, model_path = model_path_full,
    status = "ok", notes = NA_character_,
    max_rhat = diag$max_rhat, min_bulk_ess = diag$min_bulk_ess,
    min_tail_ess = diag$min_tail_ess, n_divergent = diag$n_divergent,
    max_treedepth_hits = diag$max_treedepth_hits
  )
}

grid_results <- dplyr::bind_rows(status_rows)
saveRDS(grid_results, dir_out("models", "grid_status_m07plus.rds"))
cat("\nGrid m07+ complete. Status:\n")
print(grid_results[, c("model_id", "moderator", "status", "max_rhat", "n_divergent")])
