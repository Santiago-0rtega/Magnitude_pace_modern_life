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

# Marginal means (level estimates) for location and scale.
format_location_emmeans <- function(cache, digits = 3) {
  if (is.null(cache) || is.null(cache$loc_emmeans)) return(NULL)
  .round_num(dplyr::rename_with(cache$loc_emmeans, ~ sub("^emmean$", "Estimate", .x)),
             digits)
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
