# precompute_emmeans_contrasts.R
# Run ONCE on remote-server (models are local + fast there). For each CATEGORICAL
# moderator model, compute estimated marginal means and all pairwise
# level-vs-level contrasts, for BOTH the location part and the scale (sigma)
# part, and save them as small structured tables the book can read without ever
# loading a 400 MB fit.
# `outputs/` is the compute-server workspace. Copy the completed summaries and
# tables to `Rdata/summaries/` and `Rdata/tables/` before rendering the book.
#
#   Location: emmeans(fit, ~ moderator, epred = TRUE, re_formula = NA) and
#             pairwise contrasts (no p-value adjustment). Tables and figures use
#             posterior means and equal-tail 95% credible intervals throughout.
#   Scale:    marginal means and pairwise contrasts built directly from the
#             posterior draws of the sigma (log-scale) coefficients, since
#             emmeans cannot target the scale part of a brms location-scale fit.
#             Reported on the log-sigma scale; exp() gives the residual SD.
#
# Outputs (all small):
#   outputs/summaries/emmeans_contrasts_cache.rds   (named list by model_id)
#   outputs/tables/contrasts/<id>_location_emmeans.csv
#   outputs/tables/contrasts/<id>_location_contrasts.csv
#   outputs/tables/contrasts/<id>_scale_emmeans.csv
#   outputs/tables/contrasts/<id>_scale_contrasts.csv
#
#   Rscript precompute_emmeans_contrasts.R

setwd(here::here())
source(here::here("Scripts", "01_paths.R"))
suppressMessages({
  library(brms); library(emmeans); library(dplyr); library(tibble); library(readr)
  library(coda)
})
source(here::here("Scripts", "06_model_registry.R"))

model_dir     <- dir_out("models")
tables_dir    <- dir_out("tables", "contrasts")
summ_dir      <- dir_out("summaries")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(summ_dir,   showWarnings = FALSE, recursive = TRUE)

# Prefer refit (_v2 / _v3) fits where they exist, mirroring the grid precompute.
candidates <- list(
  m06 = c("m06c_ls_taxa3", "m06b_ls_taxa_v2"),
  m09 = "m09b_ls_data_type_v2",
  m10 = "m10b_ls_transf_data_v2"
)
resolve_file <- function(mod_id, moderator) {
  cand <- c(candidates[[mod_id]], paste0(mod_id, "_ls_", moderator))
  hit  <- cand[file.exists(file.path(model_dir, paste0(cand, ".rds")))]
  if (length(hit)) hit[1] else NA_character_
}

# The actual moderator column can differ from the registry name for refit (_v2)
# models (e.g. m10 uses "transf_data_v2"). Detect it from the fit itself: it is
# the sole predictor column that is not the response or a grouping factor.
detect_moderator <- function(fit, fallback) {
  non_mod <- c("yi_lnM_safe", "ref_id", "sp_ncbi", "es_id_model", "es_id",
               "vi_lnM_safe", "es_id_db", "sys_id")
  cand <- setdiff(names(fit$data), non_mod)
  if (length(cand) == 1L) cand else fallback
}

# ── Scale (sigma) marginals + pairwise contrasts, from posterior draws ────────
# log-sigma per level = b_sigma_Intercept (+ b_sigma_<moderator><level>).
sigma_contrasts <- function(fit, moderator) {
  draws   <- brms::as_draws_df(fit)
  ref_lev <- levels(factor(fit$data[[moderator]]))[1]
  all_lev <- unique(as.character(fit$data[[moderator]]))
  all_lev <- all_lev[!is.na(all_lev)]

  sig_int <- draws[["b_sigma_Intercept"]]
  lev_draws <- list(); lev_draws[[ref_lev]] <- sig_int
  # brms strips spaces/special chars (keeps alphanumerics + underscore) when it
  # names factor-level coefficients, e.g. "Landscape change" -> "Landscapechange".
  # Match on the stripped form so levels with spaces are NOT silently dropped.
  sig_cols  <- grep(paste0("^b_sigma_", moderator), names(draws), value = TRUE)
  col_key   <- gsub("[^[:alnum:]_]", "", sub(paste0("^b_sigma_", moderator), "", sig_cols))
  for (lv in setdiff(all_lev, ref_lev)) {
    lv_key <- gsub("[^[:alnum:]_]", "", lv)
    hit    <- sig_cols[col_key == lv_key]
    if (length(hit) == 1L) {
      lev_draws[[lv]] <- sig_int + draws[[hit]]
    } else {
      message("  scale: no sigma coefficient matched level '", lv, "' — dropped")
    }
  }

  cri <- function(d) unname(quantile(d, c(0.025, 0.975)))

  emm <- bind_rows(lapply(names(lev_draws), function(lv) {
    d <- lev_draws[[lv]]; h <- cri(d)
    tibble(level = lv, emmean = mean(d), lower.CrI = h[1], upper.CrI = h[2],
           residual_SD = exp(mean(d)))
  }))

  levs  <- names(lev_draws)
  pairs <- bind_rows(lapply(seq_along(levs), function(i)
    bind_rows(lapply(seq_along(levs), function(j) {
      if (j <= i) return(NULL)
      d <- lev_draws[[levs[i]]] - lev_draws[[levs[j]]]; h <- cri(d)
      pdir <- mean(d > 0)
      tibble(contrast = paste(levs[i], "-", levs[j]),
             estimate = mean(d), lower.CrI = h[1], upper.CrI = h[2],
             pd = max(pdir, 1 - pdir))
    }))))
  list(emmeans = emm, contrasts = pairs)
}

# ── Location marginals + pairwise contrasts, via emmeans on the epred ─────────
location_contrasts <- function(fit, moderator) {
  em <- emmeans::emmeans(fit, as.formula(paste("~", moderator)),
                         epred = TRUE, re_formula = NA)
  em_draws <- as.matrix(emmeans::as.mcmc.emmGrid(em, names = FALSE))
  em_levels <- as.character(as.data.frame(em)[[1]])
  emm_df <- bind_rows(lapply(seq_along(em_levels), function(i) {
    d <- em_draws[, i]
    tibble(
      level = em_levels[i], emmean = mean(d),
      lower.CrI = unname(quantile(d, 0.025)),
      upper.CrI = unname(quantile(d, 0.975))
    )
  }))

  pw <- emmeans::contrast(em, method = "pairwise", adjust = "none")
  pw_draws <- as.matrix(emmeans::as.mcmc.emmGrid(pw, names = FALSE))
  pw_labels <- as.character(as.data.frame(pw)$contrast)
  ctr_df <- bind_rows(lapply(seq_along(pw_labels), function(i) {
    d <- pw_draws[, i]; p <- mean(d > 0)
    tibble(
      contrast = pw_labels[i], estimate = mean(d),
      lower.CrI = unname(quantile(d, 0.025)),
      upper.CrI = unname(quantile(d, 0.975)), pd = max(p, 1 - p)
    )
  }))
  list(emmeans = emm_df, contrasts = ctr_df)
}

# Explicit render-active categorical models. Do not derive this set from the
# fitting registry: m08 is an appendix result and can be absent from the primary
# grid while still requiring synchronized tables and orchard figures.
cat_grid <- tibble::tribble(
  ~model_id, ~moderator,   ~label,
  "m01",     "disturbance", "Disturbance context",
  "m02",     "design",      "Comparison design",
  "m05",     "trait_type",  "Trait type",
  "m07",     "genphen",     "Phenotypic vs genetic study",
  "m08",     "env_change",  "Environmental-change context",
  "m11",     "data_scale",  "Measurement scale"
)
cache <- list()

for (i in seq_len(nrow(cat_grid))) {
  row       <- cat_grid[i, ]
  mod_id    <- row$model_id; moderator <- row$moderator; mod_label <- row$label
  mf <- resolve_file(mod_id, moderator)
  if (is.na(mf)) { message(mod_id, ": no rds — skipped"); next }

  message(mod_id, ": reading ", mf)
  fit <- readRDS(file.path(model_dir, paste0(mf, ".rds")))
  moderator <- detect_moderator(fit, moderator)
  message("  moderator column: ", moderator)

  loc <- tryCatch(location_contrasts(fit, moderator),
                  error = function(e) { message("  location failed: ",
                                                 conditionMessage(e)); NULL })
  scl <- tryCatch(sigma_contrasts(fit, moderator),
                  error = function(e) { message("  scale failed: ",
                                                conditionMessage(e)); NULL })

  cache[[mod_id]] <- list(model_id = mod_id, moderator = moderator,
                          label = mod_label, model_file = mf,
                          loc_emmeans   = loc$emmeans,   loc_contrasts = loc$contrasts,
                          scl_emmeans   = scl$emmeans,   scl_contrasts = scl$contrasts)

  if (!is.null(loc)) {
    write_csv(loc$emmeans,   file.path(tables_dir, paste0(mod_id, "_location_emmeans.csv")))
    write_csv(loc$contrasts, file.path(tables_dir, paste0(mod_id, "_location_contrasts.csv")))
  }
  if (!is.null(scl)) {
    write_csv(scl$emmeans,   file.path(tables_dir, paste0(mod_id, "_scale_emmeans.csv")))
    write_csv(scl$contrasts, file.path(tables_dir, paste0(mod_id, "_scale_contrasts.csv")))
  }
  rm(fit); gc(verbose = FALSE)
}

attr(cache, "summary_spec") <- list(
  version = 2L, point = "posterior_mean", interval = "equal_tail_95"
)
saveRDS(cache, file.path(summ_dir, "emmeans_contrasts_cache.rds"))
message("\nDone. cached models: ", paste(names(cache), collapse = ", "))
