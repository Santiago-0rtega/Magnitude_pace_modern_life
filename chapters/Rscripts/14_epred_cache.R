# R/14_epred_cache.R
# ─────────────────────────────────────────────────────────────────────────────
# Pre-compute (once) and consume (at render) the posterior-draw caches that back
# the per-model result chapters (chapters/models/m00.qmd ... m13.qmd).
#
# WHY: loading a 400 MB brms fit and running add_epred_draws() on every render is
# slow. Instead we run the epreds ONCE with build_epred_cache() / the
# precompute_epred_draws.R runner, save a small rds per model under
# outputs/epred_draws/<id>.rds, and the chapters only read that rds and re-plot.
#
# A cache is a plain list (no brms/tidybayes classes), so chapters need only
# ggplot2 + ggdist + patchwork + dplyr to render — never brms.
#
#   $id $kind $label $moderator $predictors
#   $summary_txt   character vector  (capture.output(print(summary(fit))))
#   $fixef         matrix
#   $loc $scl      data.frames of draws (kind-dependent)
#   $raw           data.frame for continuous overlay
#   $post          data.frame for the intercept model
#   $meta          list(ndraws, n_grid, built)
#
# kind ∈ {"categorical","continuous","interaction","intercept"}
# ─────────────────────────────────────────────────────────────────────────────

## Plotting dependencies only (light — sourced by chapters at render time).
suppressMessages({
  library(ggplot2)
  library(ggdist)
  library(patchwork)
  library(dplyr)
})

# Project palette (kept in sync with R/13_epred_draws_plots.R)
COL_LOCATION <- "#4d7a5a"
COL_SCALE    <- "#e06c75"
COL_RIBBON   <- "#bcd2b1"
COL_RIBBON_S <- "#f5c2c7"
COL_POINT    <- "#675d76"

theme_cache <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      plot.title       = element_text(hjust = 0.5, face = "bold",
                                      size = rel(1.05), margin = margin(b = 6)),
      plot.subtitle    = element_text(hjust = 0.5, size = rel(0.9)),
      strip.background = element_blank(),
      legend.position  = "none"
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# BUILDERS  — run once (need brms + tidybayes; called by the runner script)
# ══════════════════════════════════════════════════════════════════════════════

# Build a newdata frame holding `keep` predictor columns and NA-ing every other
# column the model knows about, so re_formula = NA predictions never error on a
# missing grouping/id column.
.nd_fill <- function(fit, nd, keep) {
  for (col in setdiff(names(fit$data), keep)) nd[[col]] <- NA
  nd
}

build_epred_cache <- function(fit, id, kind, label,
                              moderator  = NULL,
                              predictors = NULL,
                              ndraws     = 1000,
                              n_grid     = 100) {
  stopifnot(requireNamespace("tidybayes", quietly = TRUE))
  add_epred_draws <- tidybayes::add_epred_draws

  out <- list(
    id = id, kind = kind, label = label,
    moderator = moderator, predictors = predictors,
    summary_txt = capture.output(print(summary(fit))),
    fixef = brms::fixef(fit),
    meta = list(ndraws = ndraws, n_grid = n_grid, built = Sys.time())
  )

  if (kind == "intercept") {
    dd <- brms::as_draws_df(fit)
    out$post <- data.frame(
      b_Intercept = as.numeric(dd[["b_Intercept"]]),
      sigma       = exp(as.numeric(dd[["b_sigma_Intercept"]]))
    )

  } else if (kind == "categorical") {
    nd <- fit$data |>
      dplyr::distinct(.data[[moderator]]) |>
      dplyr::filter(!is.na(.data[[moderator]]))
    nd <- .nd_fill(fit, nd, moderator)
    out$loc <- add_epred_draws(nd, fit, re_formula = NA, ndraws = ndraws) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(moderator), .draw, .epred)
    out$scl <- add_epred_draws(nd, fit, re_formula = NA, ndraws = ndraws,
                               dpar = "sigma") |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(moderator), .draw, sigma)

  } else if (kind == "continuous") {
    rng <- range(fit$data[[moderator]], na.rm = TRUE)
    nd  <- tibble::tibble(!!moderator := seq(rng[1], rng[2], length.out = n_grid))
    nd  <- .nd_fill(fit, nd, moderator)
    out$loc <- add_epred_draws(nd, fit, re_formula = NA, ndraws = ndraws) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(moderator), .row, .draw, .epred)
    out$scl <- add_epred_draws(nd, fit, re_formula = NA, ndraws = ndraws,
                               dpar = "sigma") |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(moderator), .row, .draw, sigma)
    out$raw <- fit$data[, c(moderator, "yi_lnM_safe")]

  } else if (kind == "interaction") {
    cv    <- predictors[1]
    dlevs <- levels(droplevels(factor(fit$data$disturbance)))
    rng   <- range(fit$data[[cv]], na.rm = TRUE)
    nd    <- expand.grid(x = seq(rng[1], rng[2], length.out = n_grid),
                         disturbance = dlevs, stringsAsFactors = FALSE)
    names(nd)[1] <- cv
    nd <- .nd_fill(fit, nd, c(cv, "disturbance"))
    out$loc <- add_epred_draws(nd, fit, re_formula = NA, ndraws = ndraws) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(cv), disturbance, .draw, .epred)
    out$scl <- add_epred_draws(nd, fit, re_formula = NA, ndraws = ndraws,
                               dpar = "sigma") |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(cv), disturbance, .draw, sigma)

  } else {
    stop("unknown kind: ", kind)
  }
  out
}

# ══════════════════════════════════════════════════════════════════════════════
# CONSUMERS — used by the chapters at render time (no brms / tidybayes needed)
# ══════════════════════════════════════════════════════════════════════════════

read_epred_cache <- function(id, cache_dir) {
  f <- file.path(cache_dir, paste0(id, ".rds"))
  if (file.exists(f)) readRDS(f) else NULL
}

# Print the cached summary as verbatim output, or a placeholder.
print_cache_summary <- function(cache) {
  if (is.null(cache) || is.null(cache$summary_txt)) {
    cat("Model not yet fitted — rebuild the cache once the totoro fit completes.\n")
  } else {
    writeLines(cache$summary_txt)
  }
  invisible(NULL)
}

.gg_categorical <- function(cache) {
  mod <- cache$moderator; lab <- cache$label
  p_loc <- ggplot(cache$loc, aes(x = .epred, y = .data[[mod]])) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.4) +
    stat_pointinterval(.width = c(0.50, 0.95), point_size = 2.4,
                       colour = COL_LOCATION) +
    labs(x = "Predicted lnM (location)", y = NULL, title = paste("Location:", lab)) +
    theme_cache()
  p_scl <- ggplot(cache$scl, aes(x = sigma, y = .data[[mod]])) +
    stat_pointinterval(.width = c(0.50, 0.95), point_size = 2.4,
                       colour = COL_SCALE) +
    labs(x = "Predicted sigma (residual SD)", y = NULL,
         title = paste("Scale:", lab),
         caption = "Higher sigma = more heterogeneous divergence within category") +
    theme_cache()
  (p_loc / p_scl) +
    plot_annotation(title = paste("Location–scale predictions:", lab),
                    theme = theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold", size = 13)))
}

.gg_continuous <- function(cache) {
  mod <- cache$moderator; lab <- cache$label
  p_loc <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.4)
  if (!is.null(cache$raw)) {
    p_loc <- p_loc +
      geom_point(data = cache$raw,
                 aes(x = .data[[mod]], y = yi_lnM_safe),
                 shape = 21, fill = COL_RIBBON, colour = COL_POINT,
                 size = 1.7, alpha = 0.5, stroke = 0.5)
  }
  p_loc <- p_loc +
    stat_lineribbon(data = cache$loc, aes(x = .data[[mod]], y = .epred),
                    .width = c(0.50, 0.80, 0.95), alpha = 0.30,
                    fill = COL_LOCATION, colour = NA) +
    stat_lineribbon(data = cache$loc, aes(x = .data[[mod]], y = .epred),
                    .width = 0, colour = COL_LOCATION, linewidth = 0.9) +
    labs(x = lab, y = "Predicted lnM (location)") + theme_cache()
  p_scl <- ggplot(cache$scl, aes(x = .data[[mod]], y = sigma)) +
    stat_lineribbon(.width = c(0.50, 0.80, 0.95), alpha = 0.30,
                    fill = COL_SCALE, colour = NA) +
    stat_lineribbon(.width = 0, colour = COL_SCALE, linewidth = 0.9) +
    labs(x = lab, y = "Predicted sigma (residual SD)") + theme_cache()
  (p_loc / p_scl) +
    plot_annotation(title = paste("Location–scale predictions:", lab),
                    theme = theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold", size = 13)))
}

.gg_interaction <- function(cache) {
  cv <- cache$predictors[1]; lab <- cache$label
  p_loc <- ggplot(cache$loc, aes(x = .data[[cv]], y = .epred)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.4) +
    stat_lineribbon(.width = c(0.50, 0.80, 0.95), alpha = 0.30,
                    fill = COL_LOCATION, colour = NA) +
    stat_lineribbon(.width = 0, colour = COL_LOCATION, linewidth = 0.8) +
    facet_wrap(~ disturbance) +
    labs(x = lab, y = "Predicted lnM (location)") + theme_cache()
  p_scl <- ggplot(cache$scl, aes(x = .data[[cv]], y = sigma)) +
    stat_lineribbon(.width = c(0.50, 0.80, 0.95), alpha = 0.30,
                    fill = COL_SCALE, colour = NA) +
    stat_lineribbon(.width = 0, colour = COL_SCALE, linewidth = 0.8) +
    facet_wrap(~ disturbance) +
    labs(x = lab, y = "Predicted sigma (residual SD)") + theme_cache()
  (p_loc / p_scl) +
    plot_annotation(title = paste("Location–scale predictions:", lab, "× disturbance"),
                    theme = theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold", size = 13)))
}

.gg_intercept <- function(cache) {
  p_loc <- ggplot(cache$post, aes(x = b_Intercept)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.4) +
    stat_halfeye(fill = COL_LOCATION, .width = c(0.50, 0.95), slab_alpha = 0.6) +
    labs(x = "Overall mean lnM (location)", y = NULL, title = "Grand mean") +
    theme_cache()
  p_scl <- ggplot(cache$post, aes(x = sigma)) +
    stat_halfeye(fill = COL_SCALE, .width = c(0.50, 0.95), slab_alpha = 0.6) +
    labs(x = "Baseline residual SD (exp(sigma))", y = NULL,
         title = "Residual heterogeneity") +
    theme_cache()
  (p_loc | p_scl) +
    plot_annotation(title = "Intercept-only posterior",
                    theme = theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold", size = 13)))
}

# Dispatch on cache$kind. Returns a patchwork/ggplot object, or NULL (+message).
plot_cache <- function(cache) {
  if (is.null(cache)) { cat("Figure pending model fit.\n"); return(invisible(NULL)) }
  switch(cache$kind,
    categorical = .gg_categorical(cache),
    continuous  = .gg_continuous(cache),
    interaction = .gg_interaction(cache),
    intercept   = .gg_intercept(cache),
    stop("unknown kind: ", cache$kind)
  )
}
