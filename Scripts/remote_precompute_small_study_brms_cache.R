# remote_precompute_small_study_brms_cache.R
# Run ONCE per completed fit. Builds the small per-model draw caches in
# outputs/epred_draws/small_study_<id>.rds that back 07b_small_study_effects.qmd.
# Copy completed caches to the book mirror at Rdata/epred_draws/ before render.
#
#   Rscript remote_precompute_small_study_brms_cache.R
#   FORCE=1 Rscript remote_precompute_small_study_brms_cache.R  # rebuild all
#
# Heavy step (loads each ~400 MB fit + epreds); meant to run on remote-server.
#
# Only the two models meeting the convergence criteria in
# remote_fit_one_small_study_brms.R are listed. Add n_se_sigma1 / n_v_sigmax
# here once their extended refits pass.

setwd(here::here())
source(here::here("Scripts", "01_paths.R"))
suppressMessages({
  library(brms); library(tidybayes); library(dplyr); library(tibble); library(here)
})
source(here::here("Scripts", "14_epred_cache.R"))
source(here::here("Scripts", "05_phylogeny.R"))

force_all <- nzchar(Sys.getenv("FORCE"))
model_dir <- dir_out("models")
cache_dir <- dir_out("epred_draws")
diag_dir  <- dir_out("tables", "small_study_brms")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ── Model data ───────────────────────────────────────────────────────────────
# Rebuilt with the same pipeline as remote_fit_one_small_study_brms.R so the
# bubble overlay can carry vi_lnM_safe, which the fits themselves do not store
# as a column (sampling variance enters through the es_id_model V matrix).
dat_es <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
A_full <- readRDS(dir_out("phylogeny", "proceed_A_matrix.rds"))
name_map <- readRDS(dir_out("phylogeny", "proceed_name_map.rds"))
dat_es <- apply_phylo_name_map(dat_es, name_map) |>
  dplyr::filter(
    is.finite(yi_lnM_safe), is.finite(vi_lnM_safe), vi_lnM_safe > 0,
    is.finite(n1), n1 > 0, is.finite(n2), n2 > 0
  ) |>
  dplyr::mutate(
    # Half harmonic-mean sample size: n0_tilde = n0 / 2, where the harmonic
    # mean is n0 = 2 * n1 * n2 / (n1 + n2).
    n0_tilde = (n1 * n2) / (n1 + n2),
    n_se = 1 / sqrt(n0_tilde),
    n_v = 1 / n0_tilde,
    es_id_model = factor(seq_len(dplyr::n()))
  )
dat_model <- prepare_phylo_and_data(dat_es, A_full, label = "small_study")$dat_model

# ── Builder ──────────────────────────────────────────────────────────────────
# build_epred_cache()'s continuous branch spans range(observed). The small-study
# grid must instead reach the predictor's zero bound, because that is where the
# infinite-sample, bias-adjusted intercept is read off.
.nd_fill_local <- function(fit, nd, keep) {
  for (col in setdiff(names(fit$data), keep)) nd[[col]] <- NA
  nd
}

build_small_study_cache <- function(fit, id, label, moderator, scale_label,
                                    ndraws = 1000, n_grid = 100) {
  nd <- tibble::tibble(
    !!moderator := seq(0, max(fit$data[[moderator]], na.rm = TRUE),
                       length.out = n_grid)
  )
  nd <- .nd_fill_local(fit, nd, moderator)

  loc <- tidybayes::add_epred_draws(nd, fit, re_formula = NA, ndraws = ndraws) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(moderator), .row, .draw, .epred)
  scl <- tidybayes::add_epred_draws(nd, fit, re_formula = NA, ndraws = ndraws,
                                    dpar = "sigma") |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(moderator), .row, .draw, sigma)

  diag_f <- file.path(diag_dir, paste0(id, "_diagnostics.csv"))

  list(
    id = id, kind = "continuous", label = label,
    moderator = moderator, predictors = moderator,
    scale_label = scale_label,
    summary_txt = capture.output(print(summary(fit))),
    fixef = brms::fixef(fit),
    loc = loc, scl = scl,
    raw = as.data.frame(dat_model[, c(moderator, "yi_lnM_safe", "vi_lnM_safe")]),
    diagnostics = if (file.exists(diag_f)) {
      readr::read_csv(diag_f, show_col_types = FALSE)
    } else NULL,
    meta = list(ndraws = ndraws, n_grid = n_grid, n_obs = nrow(dat_model),
                built = Sys.time())
  )
}

# ── Specs ────────────────────────────────────────────────────────────────────
# Predictors are built from the half harmonic-mean sample size,
# n0_tilde = n0 / 2 = (n1 * n2) / (n1 + n2), following Nakagawa et al. (2022).
specs <- list(
  list(id = "n_se_sigmax", moderator = "n_se",
       label = "1 / sqrt(n0_tilde)", scale_label = "sigma ~ n_se"),
  list(id = "n_v_sigma1", moderator = "n_v",
       label = "1 / n0_tilde", scale_label = "sigma ~ 1")
)

for (s in specs) {
  out_f <- file.path(cache_dir, paste0("small_study_", s$id, ".rds"))
  if (file.exists(out_f) && !force_all) {
    message(s$id, ": cache exists — skip (FORCE=1 to rebuild)"); next
  }

  fit_file <- file.path(model_dir, paste0("small_study_brms_", s$id, ".rds"))
  if (!file.exists(fit_file)) {
    message(s$id, ": no fit yet (", basename(fit_file), ") — skip"); next
  }

  message("\n", s$id, ": building from ", basename(fit_file),
          "  [", s$moderator, "]")
  fit <- readRDS(fit_file)

  cache <- build_small_study_cache(
    fit, id = s$id, label = s$label,
    moderator = s$moderator, scale_label = s$scale_label,
    ndraws = 1000, n_grid = 100
  )
  cache$summary_spec <- list(
    version = 2L, point = "posterior_mean", interval = "equal_tail_95"
  )
  cache$source_fit <- basename(fit_file)
  saveRDS(cache, out_f, compress = "xz")
  rm(fit, cache); gc(verbose = FALSE)
  message("  saved ", out_f, "  (",
          round(file.size(out_f) / 1024^2, 1), " MB)")
}

message("\nAll done.")
