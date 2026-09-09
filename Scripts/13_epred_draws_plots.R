# R/13_epred_draws_plots.R
# Bayesian posterior-draw visualizations using tidybayes + ggplot2.
# Adapted from the Collapse_environmental_predictability reference repo.
#
# Key dependencies: tidybayes, ggdist, ggplot2, patchwork, dplyr, tibble
#
# Two plot families are implemented here:
#   A) Continuous moderators  – prediction ribbons via epred_draws()
#   B) Categorical moderators – posterior dot+interval plots via epred_draws()
#
# Both families expose a location panel (mean lnM) and a scale panel (residual
# sigma), mirroring the location–scale brms model structure used throughout
# this project.
# ─────────────────────────────────────────────────────────────────────────────

library(tidybayes)
library(ggdist)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)

# ── Shared theme ──────────────────────────────────────────────────────────────

theme_draws <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.ticks.minor    = element_line(color = "black", linewidth = 0.25),
      panel.background    = element_rect(fill = "transparent", colour = NA),
      plot.background     = element_rect(fill = "transparent", colour = NA),
      axis.title.x        = element_text(margin = margin(t = 5)),
      axis.title.y        = element_text(margin = margin(r = 5)),
      plot.title          = element_text(hjust = 0.5, face = "bold",
                                         size  = rel(1.1), margin = margin(b = 6)),
      plot.subtitle       = element_text(hjust = 0.5, size = rel(0.9),
                                         margin = margin(b = 6)),
      strip.background    = element_blank(),
      legend.position     = "none"
    )
}

# Colour constants matching the project palette
COL_LOCATION <- "#4d7a5a"   # forest green  – location submodel
COL_SCALE    <- "#e06c75"   # muted red     – scale submodel
COL_RIBBON   <- "#bcd2b1"   # light green   – ribbon fill (location)
COL_RIBBON_S <- "#f5c2c7"   # light red     – ribbon fill (scale)
COL_POINT    <- "#675d76"   # muted purple  – raw data points

# ── Helper: build a newdata grid for a continuous moderator ───────────────────

.newdata_continuous <- function(fit, moderator, n_grid = 100,
                                fixed_at_zero = c("es_id_model", "ref_id", "sp_ncbi")) {
  dat  <- fit$data
  rng  <- range(dat[[moderator]], na.rm = TRUE)
  nd   <- tibble(!!moderator := seq(rng[1], rng[2], length.out = n_grid))
  for (v in fixed_at_zero) nd[[v]] <- NA
  nd
}

# ── A1) Continuous location ribbon ────────────────────────────────────────────
#
# Usage:
#   p <- plot_epred_continuous_location(
#          fit       = m03,               # fitted brms object
#          moderator = "log10_years",     # column name in fit$data
#          label     = "Elapsed time (log₁₀ years)",
#          dat_raw   = dat_model          # optional: overlay raw effect sizes
#        )

plot_epred_continuous_location <- function(fit, moderator, label,
                                           dat_raw   = NULL,
                                           ndraws    = 500,
                                           n_grid    = 100,
                                           intervals = c(0.50, 0.80, 0.95)) {
  nd <- .newdata_continuous(fit, moderator, n_grid)

  draws <- nd |>
    add_epred_draws(fit, re_formula = NA, ndraws = ndraws)

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.4)

  if (!is.null(dat_raw) && moderator %in% names(dat_raw)) {
    p <- p +
      geom_point(data = dat_raw,
                 aes(x = .data[[moderator]], y = yi_lnM_safe),
                 shape = 21, fill = COL_RIBBON, colour = COL_POINT,
                 size = 1.8, alpha = 0.55, stroke = 0.6)
  }

  p +
    stat_lineribbon(data = draws,
                    aes(x = .data[[moderator]], y = .epred),
                    .width = intervals,
                    alpha  = 0.30,
                    fill   = COL_LOCATION,
                    colour = NA) +
    stat_lineribbon(data = draws,
                    aes(x = .data[[moderator]], y = .epred),
                    .width = 0,              # median line only
                    colour = COL_LOCATION,
                    linewidth = 0.9,
                    fill   = NA) +
    labs(
      x        = label,
      y        = "Predicted lnM (location)",
      subtitle = paste0(intervals * 100, "% CrI", collapse = ", ")
    ) +
    theme_draws()
}

# ── A2) Continuous scale ribbon ───────────────────────────────────────────────
#
# Extracts the sigma distributional parameter from a location–scale brms model.

plot_epred_continuous_scale <- function(fit, moderator, label,
                                        ndraws    = 500,
                                        n_grid    = 100,
                                        intervals = c(0.50, 0.80, 0.95)) {
  nd <- .newdata_continuous(fit, moderator, n_grid)

  draws <- nd |>
    add_epred_draws(fit, re_formula = NA, ndraws = ndraws, dpar = "sigma")

  ggplot() +
    stat_lineribbon(data = draws,
                    aes(x = .data[[moderator]], y = sigma),
                    .width = intervals,
                    alpha  = 0.30,
                    fill   = COL_SCALE,
                    colour = NA) +
    stat_lineribbon(data = draws,
                    aes(x = .data[[moderator]], y = sigma),
                    .width    = 0,
                    colour    = COL_SCALE,
                    linewidth = 0.9,
                    fill      = NA) +
    labs(
      x        = label,
      y        = "Predicted σ (residual SD)",
      caption  = "Higher σ = greater residual heterogeneity in divergence"
    ) +
    theme_draws()
}

# ── A3) Combined continuous location + scale (patchwork) ─────────────────────

plot_epred_continuous_ls <- function(fit, moderator, label,
                                     dat_raw = NULL, ...) {
  p_loc <- plot_epred_continuous_location(fit, moderator, label,
                                          dat_raw = dat_raw, ...)
  p_scl <- plot_epred_continuous_scale(fit, moderator, label, ...)

  (p_loc / p_scl) +
    plot_annotation(
      title = paste("Location–scale predictions:", label),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 13))
    )
}

# ── B1) Categorical location – posterior draws as dot+interval ────────────────
#
# Uses epred_draws() to show the full posterior for each category level.
#
# Usage:
#   p <- plot_epred_categorical_location(
#          fit       = m01,
#          moderator = "disturbance",
#          label     = "Disturbance context"
#        )

plot_epred_categorical_location <- function(fit, moderator, label,
                                             ndraws    = 1000,
                                             intervals = c(0.50, 0.95)) {
  nd <- fit$data |>
    distinct(.data[[moderator]]) |>
    filter(!is.na(.data[[moderator]])) |>
    mutate(es_id_model = NA, ref_id = NA, sp_ncbi = NA)

  draws <- nd |>
    add_epred_draws(fit, re_formula = NA, ndraws = ndraws)

  ggplot(draws, aes(x = .epred, y = .data[[moderator]])) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.4) +
    stat_pointinterval(
      aes(colour = after_stat(level)),
      point_size   = 2.5,
      .width       = intervals,
      point_colour = COL_LOCATION
    ) +
    scale_colour_manual(
      values = setNames(
        colorRampPalette(c(COL_RIBBON, COL_LOCATION))(length(intervals)),
        as.character(intervals)
      ),
      guide = "none"
    ) +
    labs(
      x     = "Predicted lnM (location)",
      y     = NULL,
      title = paste("Location:", label)
    ) +
    theme_draws()
}

# ── B2) Categorical scale – posterior draws as dot+interval ───────────────────

plot_epred_categorical_scale <- function(fit, moderator, label,
                                          ndraws    = 1000,
                                          intervals = c(0.50, 0.95)) {
  nd <- fit$data |>
    distinct(.data[[moderator]]) |>
    filter(!is.na(.data[[moderator]])) |>
    mutate(es_id_model = NA, ref_id = NA, sp_ncbi = NA)

  draws <- nd |>
    add_epred_draws(fit, re_formula = NA, ndraws = ndraws, dpar = "sigma")

  ggplot(draws, aes(x = sigma, y = .data[[moderator]])) +
    stat_pointinterval(
      aes(colour = after_stat(level)),
      point_size   = 2.5,
      .width       = intervals,
      point_colour = COL_SCALE
    ) +
    scale_colour_manual(
      values = setNames(
        colorRampPalette(c(COL_RIBBON_S, COL_SCALE))(length(intervals)),
        as.character(intervals)
      ),
      guide = "none"
    ) +
    labs(
      x       = "Predicted σ (residual SD)",
      y       = NULL,
      title   = paste("Scale:", label),
      caption = "Higher σ = more heterogeneous divergence within category"
    ) +
    theme_draws()
}

# ── B3) Combined categorical location + scale ─────────────────────────────────

plot_epred_categorical_ls <- function(fit, moderator, label, ...) {
  p_loc <- plot_epred_categorical_location(fit, moderator, label, ...)
  p_scl <- plot_epred_categorical_scale(fit, moderator, label, ...)

  (p_loc / p_scl) +
    plot_annotation(
      title = paste("Location–scale predictions:", label),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 13))
    )
}

# ── C) Full grid: iterate over moderator_grid ─────────────────────────────────
#
# Convenience wrapper that loops over the moderator_grid tibble and calls the
# right plot family for each model, then saves via save_plot_dual().
#
# Requires: moderator_grid, fit_or_read_model(), save_plot_dual()
#           (sourced from R/06_model_registry.R, R/08_fit_or_read_model.R,
#                        R/12_save_figures.R)
#
# Usage (from a .qmd chunk or script):
#   source(here::here("Scripts", "13_epred_draws_plots.R"))
#   generate_all_epred_plots(dat_model = proceed_lnm)

generate_all_epred_plots <- function(dat_model,
                                     ndraws    = 500,
                                     intervals = c(0.50, 0.80, 0.95),
                                     refit     = FALSE) {
  stopifnot(exists("moderator_grid"), exists("fit_or_read_model"),
            exists("save_plot_dual"))

  for (i in seq_len(nrow(moderator_grid))) {
    row  <- moderator_grid[i, ]
    mid  <- row$model_id
    mod  <- row$moderator
    lbl  <- row$label
    type <- row$type

    message("epred plot: ", mid, " [", type, "] — ", lbl)

    fit <- tryCatch(
      fit_or_read_model(mid, refit = refit),
      error = function(e) {
        message("  skipped (model not found): ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(fit)) next

    if (type == "continuous") {
      p <- plot_epred_continuous_ls(fit, mod, lbl,
                                    dat_raw   = dat_model,
                                    ndraws    = ndraws,
                                    intervals = intervals)
    } else {
      p <- plot_epred_categorical_ls(fit, mod, lbl,
                                     ndraws    = ndraws,
                                     intervals = intervals)
    }

    save_plot_dual(p, paste0(mid, "_", mod, "_epred_ls"),
                   width = 180, height = 260)
    message("  saved: ", mid, "_", mod, "_epred_ls")
  }

  invisible(NULL)
}
