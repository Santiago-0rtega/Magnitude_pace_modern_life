# precompute_grid_summaries.R
# Run ONCE on totoro (models are local + fast there). Loads each moderator-grid
# model, extracts fixed effects + diagnostics + status, and writes the small
# tables that 05_location_scale_model_grid.qmd and 08_publication_figures.qmd
# consume — so the BOOK never loads a 400 MB model.
#
# Outputs (all small) under outputs/tables/:
#   location_effects_all_moderators.csv
#   scale_effects_all_moderators.csv
#   model_diagnostics_all_moderators.csv
#   model_grid_status.csv
#   grid_summary_cache.rds   (bundle: list(status, location, scale, diagnostics))
#
#   Rscript precompute_grid_summaries.R

setwd(if (dir.exists("/home/ortegara/Documents/PACE"))
        "/home/ortegara/Documents/PACE" else getwd())
suppressMessages({ library(brms); library(dplyr); library(tibble); library(readr); library(purrr) })
source(file.path("R", "00_packages.R"))
source(file.path("R", "01_paths.R"))
source(file.path("R", "06_model_registry.R"))
source(file.path("R", "09_model_summaries.R"))
source(file.path("R", "10_model_diagnostics.R"))

model_dir  <- file.path("outputs", "models")
tables_dir <- file.path("outputs", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

dat_es <- readRDS(file.path("outputs", "effect_sizes", "proceed_lnm_safe.rds"))

# Candidate rds per model_id (prefer the *_v2 refit, fall back to base) —
# mirrors the per-model chapter / epred-cache resolution.
# m06/m09 have no entry here — they're excluded from moderator_grid (see
# R/06_model_registry.R) so this loop never resolves a file for them.
candidates <- list(
  m10 = "m10b_ls_transf_data_v2"
)
resolve_file <- function(mod_id, moderator) {
  cand <- c(candidates[[mod_id]], paste0(mod_id, "_ls_", moderator))
  hit  <- cand[file.exists(file.path(model_dir, paste0(cand, ".rds")))]
  if (length(hit)) hit[1] else NA_character_
}

all_fe <- list(); all_diag <- list(); status_rows <- list()

for (i in seq_len(nrow(moderator_grid))) {
  row       <- moderator_grid[i, ]
  mod_id    <- row$model_id; moderator <- row$moderator
  mod_label <- row$label;    mod_type  <- row$type

  na_row <- function(status, notes, nr = 0L, nl = NA_integer_, nx = NA_integer_,
                     path = NA_character_)
    tibble(model_id = mod_id, moderator = moderator, label = mod_label,
           type = mod_type, n_rows = nr, n_levels = nl, n_missing_excluded = nx,
           model_path = path, status = status, notes = notes,
           max_rhat = NA_real_, min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
           n_divergent = NA_integer_, max_treedepth_hits = NA_integer_)

  if (!moderator %in% names(dat_es)) {
    status_rows[[i]] <- na_row("skipped: variable absent", "Moderator not found.")
    next
  }
  dat_model <- dat_es[!is.na(dat_es[[moderator]]), ]

  mf <- resolve_file(mod_id, moderator)
  if (is.na(mf)) {
    status_rows[[i]] <- na_row("missing: no cached fit",
                               "No model .rds found for this moderator.",
                               nr = nrow(dat_model), nl = n_levels, nx = n_missing)
    message(mod_id, ": no rds — skipped"); next
  }

  message(mod_id, ": reading ", mf)
  fit  <- readRDS(file.path(model_dir, paste0(mf, ".rds")))
  fit_data <- fit$data
  n_missing <- nrow(dat_es) - nrow(fit_data)
  n_levels  <- if (mod_type == "categorical")
    length(unique(droplevels(factor(fit_data[[moderator]])))) else NA_integer_
  fe   <- extract_fixed_effects(fit, mod_id, moderator)
  diag <- tryCatch(
    extract_diagnostics(fit, mod_id, moderator),
    error = function(e) tibble(
      model_id = mod_id, moderator = moderator,
      max_rhat = NA_real_, min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
      n_divergent = NA_integer_, max_treedepth_hits = NA_integer_))

  all_fe[[mod_id]]   <- split_location_scale(fe)
  all_diag[[mod_id]] <- diag
  status_rows[[i]] <- tibble(
    model_id = mod_id, moderator = moderator, label = mod_label, type = mod_type,
    n_rows = nrow(fit_data), n_levels = n_levels, n_missing_excluded = n_missing,
    model_path = file.path(model_dir, paste0(mf, ".rds")), status = "fitted",
    notes = if (startsWith(mf, paste0(mod_id, "b"))) "refit (_v2)" else "",
    max_rhat = diag$max_rhat, min_bulk_ess = diag$min_bulk_ess,
    min_tail_ess = diag$min_tail_ess, n_divergent = diag$n_divergent,
    max_treedepth_hits = diag$max_treedepth_hits)
  rm(fit); gc(verbose = FALSE)
}

model_status <- bind_rows(status_rows)
summ <- save_combined_summaries(all_fe, tables_dir = tables_dir)
dtbl <- save_combined_diagnostics(all_diag, tables_dir = tables_dir)
write_csv(model_status, file.path(tables_dir, "model_grid_status.csv"))

saveRDS(list(status = model_status, location = summ$location,
             scale = summ$scale, diagnostics = dtbl),
        file.path(tables_dir, "grid_summary_cache.rds"))

message("\nDone. fitted models: ", sum(model_status$status == "fitted"),
        " / ", nrow(model_status))
