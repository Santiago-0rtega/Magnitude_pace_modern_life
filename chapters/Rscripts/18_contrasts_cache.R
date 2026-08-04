# R/18_contrasts_cache.R
# ─────────────────────────────────────────────────────────────────────────────
# Read-side helpers for the emmeans contrast tables and the m00 heterogeneity
# decomposition. Chapters consume the small caches produced on totoro by
# precompute_emmeans_contrasts.R and precompute_heterogeneity.R — no brms /
# emmeans / 400 MB fit is ever loaded at render time.
#
#   outputs/summaries/emmeans_contrasts_cache.rds   list by model_id, each:
#     $loc_emmeans $loc_contrasts $scl_emmeans $scl_contrasts
#   outputs/summaries/heterogeneity_m00.rds         list($I2, $components, $Vbar)
# ─────────────────────────────────────────────────────────────────────────────

suppressMessages({ library(dplyr); library(kableExtra) })

# CI excludes 0 (both bounds same sign) → contrast is "credibly non-zero".
.ci_excludes_zero <- function(lo, hi) (lo > 0 & hi > 0) | (lo < 0 & hi < 0)

read_contrasts_cache <- function(
    id,
    path = here::here("Rdata", "summaries", "emmeans_contrasts_cache.rds")) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)[[id]]
}

format_heterogeneity_components <- function(het) {
  if (is.null(het) || is.null(het$components)) return(NULL)
  df <- het$components
  tibble::tibble(
    Quantity = df$quantity,
    Value    = df$value
  )
}

# Round only the numeric columns, leave labels alone.
.round_num <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return(df)
  dplyr::mutate(df, dplyr::across(dplyr::where(is.numeric), ~ round(.x, digits)))
}

# Location pairwise contrasts, tidied for display.
# emmeans summary columns are typically: contrast, estimate, lower.HPD,
# upper.HPD (Bayesian) — plus our added pd. We rename to friendly headers.
format_location_contrasts <- function(cache, digits = 3) {
  if (is.null(cache) || is.null(cache$loc_contrasts)) return(NULL)
  df <- cache$loc_contrasts
  lo <- dplyr::coalesce(df[["lower.HPD"]], df[["asymp.LCL"]], df[["lower.CL"]])
  hi <- dplyr::coalesce(df[["upper.HPD"]], df[["asymp.UCL"]], df[["upper.CL"]])
  out <- tibble::tibble(
    Contrast          = as.character(df$contrast),
    Estimate          = df$estimate,
    `Lower 95% HPD`   = lo,
    `Upper 95% HPD`   = hi
  )
  out[[".sig"]] <- .ci_excludes_zero(lo, hi)
  .round_num(out, digits)
}

# Scale (sigma) pairwise contrasts, on the log-sigma scale.
format_scale_contrasts <- function(cache, digits = 3) {
  if (is.null(cache) || is.null(cache$scl_contrasts)) return(NULL)
  df <- cache$scl_contrasts
  out <- tibble::tibble(
    Contrast          = as.character(df$contrast),
    `Estimate (logσ)` = df$estimate,
    `Lower 95% HPD`   = df$lower.HPD,
    `Upper 95% HPD`   = df$upper.HPD
  )
  out[[".sig"]] <- .ci_excludes_zero(df$lower.HPD, df$upper.HPD)
  .round_num(out, digits)
}

# Convert a location estimate on the lnM scale to the approximate equivalent
# standardized mean difference described in the lnM chapter. The transform is
# strictly increasing, so applying it to both CrI endpoints preserves coverage.
lnm_to_d_eq <- function(lnm) sqrt(2 * exp(2 * lnm))

# Summarize posterior draws on the lnM and equivalent-d scales. Quantiles are
# calculated from the draws; because the transform is monotonic, transformed
# endpoint quantiles are identical to quantiles of the transformed draws.
format_posterior_d_eq <- function(lnm_draws, label = "Overall", digits = 3) {
  if (is.null(lnm_draws) || !length(lnm_draws)) return(NULL)
  lnm_est <- mean(lnm_draws, na.rm = TRUE)
  lnm_ci  <- stats::quantile(lnm_draws, c(0.025, 0.975), na.rm = TRUE)
  .round_num(tibble::tibble(
    Estimate              = label,
    `Mean lnM`            = lnm_est,
    `lnM lower 95% CrI`   = unname(lnm_ci[1]),
    `lnM upper 95% CrI`   = unname(lnm_ci[2]),
    `d_eq`                = lnm_to_d_eq(lnm_est),
    `d_eq lower 95% CrI`  = lnm_to_d_eq(unname(lnm_ci[1])),
    `d_eq upper 95% CrI`  = lnm_to_d_eq(unname(lnm_ci[2]))
  ), digits)
}

# Summarize cached continuous-model predictions at interpretable values on the
# original time scale (1, 10, and 100 years or generations by default).
format_continuous_d_eq <- function(cache, original_values = c(1, 10, 100),
                                   unit = "Time", digits = 3) {
  if (is.null(cache) || is.null(cache$loc) || is.null(cache$moderator)) return(NULL)
  mod <- cache$moderator
  target_log10 <- log10(original_values)
  available <- sort(unique(cache$loc[[mod]]))
  keep <- target_log10 >= min(available) & target_log10 <= max(available)
  target_log10 <- target_log10[keep]
  original_values <- original_values[keep]
  if (!length(target_log10)) return(NULL)

  rows <- lapply(seq_along(target_log10), function(i) {
    grid_value <- available[which.min(abs(available - target_log10[i]))]
    draws <- cache$loc$.epred[cache$loc[[mod]] == grid_value]
    out <- format_posterior_d_eq(draws, label = as.character(original_values[i]),
                                 digits = digits)
    names(out)[1] <- unit
    dplyr::mutate(out, `log10 value` = grid_value, .after = 1)
  })
  dplyr::bind_rows(rows)
}

# Marginal means (level estimates) for location and scale.
format_location_emmeans <- function(cache, digits = 3) {
  if (is.null(cache) || is.null(cache$loc_emmeans)) return(NULL)
  df <- cache$loc_emmeans
  lo_name <- intersect(c("lower.HPD", "asymp.LCL", "lower.CL"), names(df))[1]
  hi_name <- intersect(c("upper.HPD", "asymp.UCL", "upper.CL"), names(df))[1]
  if (is.na(lo_name) || is.na(hi_name)) return(NULL)
  lo <- df[[lo_name]]
  hi <- df[[hi_name]]

  out <- tibble::tibble(
    Level                = as.character(df$level),
    `Mean lnM`           = df$emmean,
    `lnM lower 95% CrI`  = lo,
    `lnM upper 95% CrI`  = hi,
    `d_eq`               = lnm_to_d_eq(df$emmean),
    `d_eq lower 95% CrI` = lnm_to_d_eq(lo),
    `d_eq upper 95% CrI` = lnm_to_d_eq(hi)
  )
  .round_num(out, digits)
}
format_scale_emmeans <- function(cache, digits = 3) {
  if (is.null(cache) || is.null(cache$scl_emmeans)) return(NULL)
  .round_num(cache$scl_emmeans, digits)
}

# knitr::kable wrapper that prints a placeholder instead of erroring on NULL.
kable_contrasts <- function(df, caption = NULL, digits = 3) {
  if (is.null(df) || !nrow(df))
    return(knitr::asis_output("_Contrast table unavailable — rebuild the cache on totoro._"))
  sig <- if (".sig" %in% names(df)) df[[".sig"]] else rep(FALSE, nrow(df))
  df  <- df[, setdiff(names(df), ".sig"), drop = FALSE]     # hide the flag column
  kb  <- kableExtra::kbl(df, caption = caption, digits = digits) |>
    kableExtra::kable_styling(full_width = FALSE)
  rows <- which(sig)                                        # NA treated as FALSE
  if (length(rows)) kb <- kableExtra::row_spec(kb, rows, bold = TRUE)
  kb
}

# ── Heterogeneity (m00) ───────────────────────────────────────────────────────
read_heterogeneity <- function(
    path = here::here("Rdata", "summaries", "heterogeneity_m00.rds")) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

format_heterogeneity <- function(het) {
  if (is.null(het) || is.null(het$I2)) return(NULL)
  df <- het$I2
  tibble::tibble(
    Component            = df$component,
    `I² (%)`             = df$median,
    `Lower 95% CrI`      = df$l95,
    `Upper 95% CrI`      = df$u95
  )
}
