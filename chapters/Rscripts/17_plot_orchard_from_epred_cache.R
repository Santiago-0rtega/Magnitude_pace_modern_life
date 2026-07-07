library(here)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)

out_dir <- here::here("Rdata", "figures", "publication", "orchard")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat_es <- readRDS(here::here("Rdata", "effect_sizes", "proceed_lnm_safe.rds"))
cache_dir <- here::here("Rdata", "epred_draws")

model_ids <- Sys.getenv("ORCHARD_MODEL_IDS", unset = "")
model_ids <- if (nzchar(model_ids)) {
  stringr::str_split(model_ids, "\\s*,\\s*")[[1]]
} else {
  sprintf("m%02d", 0:11)
}

label_maps <- list(
  m01 = c(
    "Climate change" = "Climate change",
    "Hunt_harv" = "Hunting / harvest",
    "Introduction" = "Introduction",
    "Landscape change" = "Landscape change",
    "Other" = "Other (in situ natural variation)",
    "Pollution" = "Pollution",
    "Response to introductions" = "Response to introductions"
  ),
  m02 = c("Allochronic" = "Allochronic", "Synchronic" = "Synchronic"),
  m05 = c(
    "behavior" = "Behavior", "growth" = "Growth", "lifehistory" = "Life history",
    "otherLH" = "Other life history", "othermorphology" = "Other morphology",
    "phenology" = "Phenology", "physio" = "Physiology", "size" = "Size"
  ),
  m06 = c(
    "Arthropod" = "Arthropod", "Bird" = "Bird", "Fish" = "Fish",
    "Mammal" = "Mammal", "Mollusk" = "Mollusk", "Plant" = "Plant",
    "Reptile" = "Reptile"
  ),
  m07 = c("Genetic" = "Genetic", "Phenotypic" = "Phenotypic"),
  m08 = c("novel" = "Novel", "ongoing" = "Ongoing"),
  m09 = c(
    "area" = "Area", "count" = "Count", "cube (3D)" = "Cube (3D)",
    "linear" = "Linear", "mass" = "Mass", "proportion" = "Proportion",
    "rate" = "Rate", "time" = "Time", "volume" = "Volume"
  ),
  m10 = c(
    "raw" = "Raw", "ln" = "ln", "log10" = "log10", "arcsin" = "arcsin",
    "resid" = "Residuals", "ord" = "Ordinal", "other" = "Other transformed",
    "Transformed" = "Transformed (ln / log10 / arcsin / resid / ord)"
  ),
  m11 = c("interval" = "Interval (arbitrary zero)", "ratio" = "Ratio (true zero)")
)

COL_LOCATION <- "#0072B2"
COL_LOCATION_LIGHT <- "#88CCEE"
COL_SCALE <- "#D55E00"
COL_SCALE_LIGHT <- "#E69F00"
COL_SCALE_RIBBON <- "#F2B27E"
PRECISION_LABEL <- "Effect-size precision (1/SE)"
SCALE_CAPTION_MARKS <- "Points show absolute residual lnM values; diamonds and intervals show model-estimated residual heterogeneity (sigma_lnm)."
SCALE_CAPTION_RIBBON <- "Points show absolute residual lnM values; line and ribbon show model-estimated residual heterogeneity (sigma_lnm)."
CATEGORY_COLS <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
  "#56B4E9", "#6A3D9A", "#999999", "#000000"
)

theme_orchard <- function() {
  ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_text(size = 11),
      axis.title.x       = ggplot2::element_text(size = 12),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey92"),
      legend.position    = "bottom",
      legend.title       = ggplot2::element_text(size = 10),
      plot.title         = ggplot2::element_text(size = 13, face = "bold"),
      plot.caption       = ggplot2::element_text(size = 9, colour = "grey50")
    )
}

save_plot <- function(p, stem, width = 9, height = 6) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
                  width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), p,
                  width = width, height = height, dpi = 300, type = "cairo")
}

apply_labels <- function(x, map) {
  x <- as.character(x)
  if (is.null(map)) return(x)
  dplyr::recode(x, !!!map, .default = x)
}

summarise_draws <- function(data, group, value) {
  data |>
    dplyr::group_by(.data[[group]]) |>
    dplyr::summarise(
      estimate = median(.data[[value]], na.rm = TRUE),
      lowerCL = quantile(.data[[value]], 0.025, na.rm = TRUE),
      upperCL = quantile(.data[[value]], 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_intercept <- function(cache) {
  raw <- dat_es |> dplyr::mutate(precision = 1 / sqrt(vi_lnM_safe))
  int_est <- median(cache$post$b_Intercept)
  int_lo <- quantile(cache$post$b_Intercept, 0.025)
  int_hi <- quantile(cache$post$b_Intercept, 0.975)
  sig_est <- median(cache$post$sigma)
  sig_lo <- quantile(cache$post$sigma, 0.025)
  sig_hi <- quantile(cache$post$sigma, 0.975)
  cap_half_height <- 0.08

  p_loc <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggbeeswarm::geom_quasirandom(
      data = raw,
      ggplot2::aes(x = yi_lnM_safe, y = 1, size = precision),
      alpha = 0.20, shape = 21, fill = COL_LOCATION_LIGHT, colour = COL_LOCATION,
      groupOnX = FALSE
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = int_lo, xend = int_hi, y = 1, yend = 1),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = int_lo, xend = int_lo,
                   y = 1 - cap_half_height, yend = 1 + cap_half_height),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = int_hi, xend = int_hi,
                   y = 1 - cap_half_height, yend = 1 + cap_half_height),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_point(ggplot2::aes(x = int_est, y = 1),
                        size = 3.4, shape = 23,
                        fill = "white", colour = "grey10", stroke = 1.1) +
    ggplot2::scale_size_continuous(name = PRECISION_LABEL, range = c(0.3, 4)) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::labs(x = "Location effect (lnM)", y = NULL,
                  title = "Location -- Overall baseline (m00)") +
    theme_orchard()

  raw_sig <- raw |> dplyr::mutate(abs_resid = abs(yi_lnM_safe - int_est))
  p_scl <- ggplot2::ggplot() +
    ggbeeswarm::geom_quasirandom(
      data = raw_sig,
      ggplot2::aes(x = abs_resid, y = 1, size = precision),
      alpha = 0.20, shape = 21, fill = COL_SCALE_LIGHT, colour = COL_SCALE,
      groupOnX = FALSE
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = sig_lo, xend = sig_hi, y = 1, yend = 1),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = sig_lo, xend = sig_lo,
                   y = 1 - cap_half_height, yend = 1 + cap_half_height),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = sig_hi, xend = sig_hi,
                   y = 1 - cap_half_height, yend = 1 + cap_half_height),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_point(ggplot2::aes(x = sig_est, y = 1),
                        size = 3.4, shape = 23,
                        fill = "white", colour = "grey10", stroke = 1.1) +
    ggplot2::scale_size_continuous(name = PRECISION_LABEL, range = c(0.3, 4)) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::labs(x = "residual lnM (SD)", y = NULL,
                  title = "Scale -- Overall baseline (m00)",
                  caption = SCALE_CAPTION_MARKS) +
    theme_orchard()

  list(location = p_loc, scale = p_scl,
       combined = patchwork::wrap_plots(p_loc, p_scl, ncol = 1, guides = "collect") &
         ggplot2::theme(legend.position = "bottom"),
       width = 9, height = 4, combined_height = 8)
}

plot_continuous <- function(cache) {
  mod <- cache$moderator
  lab <- cache$label
  raw <- dat_es[!is.na(dat_es[[mod]]), ] |>
    dplyr::mutate(precision = 1 / sqrt(vi_lnM_safe))
  pred <- summarise_draws(cache$loc, mod, ".epred")
  pred_sig <- summarise_draws(cache$scl, mod, "sigma")
  raw_pred_loc <- approx(pred[[mod]], pred$estimate, xout = raw[[mod]], rule = 2)$y
  raw_sig <- raw |> dplyr::mutate(abs_resid = abs(yi_lnM_safe - raw_pred_loc))

  p_loc <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(
      data = raw,
      ggplot2::aes(x = .data[[mod]], y = yi_lnM_safe, size = precision),
      alpha = 0.22, shape = 21, fill = COL_LOCATION_LIGHT, colour = COL_LOCATION
    ) +
    ggplot2::geom_ribbon(
      data = pred,
      ggplot2::aes(x = .data[[mod]], ymin = lowerCL, ymax = upperCL),
      alpha = 0.35, fill = COL_LOCATION
    ) +
    ggplot2::geom_line(
      data = pred,
      ggplot2::aes(x = .data[[mod]], y = estimate),
      linewidth = 1.1, colour = COL_LOCATION
    ) +
    ggplot2::scale_size_continuous(name = PRECISION_LABEL, range = c(0.3, 4)) +
    ggplot2::labs(x = lab, y = "lnM", title = paste("Location --", lab)) +
    theme_orchard()

  p_scl <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = raw_sig,
      ggplot2::aes(x = .data[[mod]], y = abs_resid, size = precision),
      alpha = 0.20, shape = 21, fill = COL_SCALE_LIGHT, colour = COL_SCALE
    ) +
    ggplot2::geom_ribbon(
      data = pred_sig,
      ggplot2::aes(x = .data[[mod]], ymin = lowerCL, ymax = upperCL),
      alpha = 0.35, fill = COL_SCALE_RIBBON
    ) +
    ggplot2::geom_line(
      data = pred_sig,
      ggplot2::aes(x = .data[[mod]], y = estimate),
      linewidth = 1.1, colour = COL_SCALE
    ) +
    ggplot2::scale_size_continuous(name = PRECISION_LABEL, range = c(0.3, 4)) +
    ggplot2::labs(x = lab, y = "residual lnM (SD)",
                  title = paste("Scale --", lab),
                  caption = SCALE_CAPTION_RIBBON) +
    theme_orchard()

  list(location = p_loc, scale = p_scl,
       combined = patchwork::wrap_plots(p_loc, p_scl, ncol = 1, guides = "collect") &
         ggplot2::theme(legend.position = "bottom"),
       width = 8, height = 5, combined_height = 9)
}

plot_categorical <- function(cache) {
  mod <- cache$moderator
  lab <- cache$label
  map <- label_maps[[cache$id]]

  ests <- summarise_draws(cache$loc, mod, ".epred") |>
    dplyr::mutate(level = apply_labels(.data[[mod]], map))
  sig_ests <- summarise_draws(cache$scl, mod, "sigma") |>
    dplyr::mutate(level = apply_labels(.data[[mod]], map))

  raw <- if (mod %in% names(dat_es)) {
    dat_es[!is.na(dat_es[[mod]]), ] |>
      dplyr::mutate(
        level = apply_labels(.data[[mod]], map),
        precision = 1 / sqrt(vi_lnM_safe)
      )
  } else {
    NULL
  }

  lev_order <- ests$level[order(ests$estimate)]
  ests$level <- factor(ests$level, levels = lev_order)
  sig_ests$level <- factor(sig_ests$level, levels = lev_order)
  if (!is.null(raw)) raw$level <- factor(raw$level, levels = lev_order)
  n_levels <- length(lev_order)
  level_cols <- CATEGORY_COLS[seq_len(n_levels)]
  names(level_cols) <- lev_order
  y_breaks <- seq_len(n_levels)
  ests <- ests |> dplyr::mutate(y_index = as.numeric(level))
  sig_ests <- sig_ests |> dplyr::mutate(y_index = as.numeric(level))
  if (!is.null(raw)) raw <- raw |> dplyr::mutate(y_index = as.numeric(level))
  cap_half_height <- 0.12

  p_loc <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    { if (!is.null(raw)) ggbeeswarm::geom_quasirandom(
      data = raw,
      ggplot2::aes(x = yi_lnM_safe, y = y_index,
                   size = precision, colour = level, fill = level),
      alpha = 0.30, shape = 21, groupOnX = FALSE
    ) } +
    ggplot2::geom_segment(
      data = ests,
      ggplot2::aes(x = lowerCL, xend = upperCL, y = y_index, yend = y_index),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_segment(
      data = ests,
      ggplot2::aes(x = lowerCL, xend = lowerCL,
                   y = y_index - cap_half_height, yend = y_index + cap_half_height),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_segment(
      data = ests,
      ggplot2::aes(x = upperCL, xend = upperCL,
                   y = y_index - cap_half_height, yend = y_index + cap_half_height),
      linewidth = 1.45, colour = "grey10"
    ) +
    ggplot2::geom_point(
      data = ests,
      ggplot2::aes(x = estimate, y = y_index),
      size = 3.4, shape = 23, fill = "white", colour = "grey10", stroke = 1.1
    ) +
    ggplot2::scale_colour_manual(values = level_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = level_cols, guide = "none") +
    ggplot2::scale_size_continuous(name = PRECISION_LABEL, range = c(0.4, 4)) +
    ggplot2::scale_y_continuous(breaks = y_breaks, labels = lev_order) +
    ggplot2::labs(x = "Location effect (lnM)", y = NULL,
                  title = paste("Location --", lab)) +
    theme_orchard()

  raw_sig <- if (!is.null(raw)) {
    raw |>
      dplyr::left_join(dplyr::select(ests, level, loc_est = estimate), by = "level") |>
      dplyr::mutate(abs_resid = abs(yi_lnM_safe - loc_est))
  } else {
    NULL
  }

  p_scl <- ggplot2::ggplot() +
    { if (!is.null(raw_sig)) ggbeeswarm::geom_quasirandom(
      data = raw_sig,
      ggplot2::aes(x = abs_resid, y = y_index,
                   size = precision, colour = level, fill = level),
      alpha = 0.30, shape = 21, groupOnX = FALSE
    ) } +
    ggplot2::geom_segment(
      data = sig_ests,
      ggplot2::aes(x = lowerCL, xend = upperCL, y = y_index, yend = y_index),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_segment(
      data = sig_ests,
      ggplot2::aes(x = lowerCL, xend = lowerCL,
                   y = y_index - cap_half_height, yend = y_index + cap_half_height),
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_segment(
      data = sig_ests,
      ggplot2::aes(x = upperCL, xend = upperCL,
                   y = y_index - cap_half_height, yend = y_index + cap_half_height),
      linewidth = 1.45, colour = "grey10"
    ) +
    ggplot2::geom_point(
      data = sig_ests,
      ggplot2::aes(x = estimate, y = y_index),
      size = 3.4, shape = 23, fill = "white", colour = "grey10", stroke = 1.1
    ) +
    ggplot2::scale_colour_manual(values = level_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = level_cols, guide = "none") +
    ggplot2::scale_size_continuous(name = PRECISION_LABEL, range = c(0.4, 4)) +
    ggplot2::scale_y_continuous(breaks = y_breaks, labels = lev_order) +
    ggplot2::labs(x = "residual lnM (SD)", y = NULL,
                  title = paste("Scale --", lab),
                  caption = SCALE_CAPTION_MARKS) +
    theme_orchard()

  height_single <- max(2.5 + n_levels * 0.55, 5)
  list(location = p_loc, scale = p_scl,
       combined = patchwork::wrap_plots(p_loc, p_scl, ncol = 1, guides = "collect") &
         ggplot2::theme(legend.position = "bottom"),
       width = 9, height = height_single,
       combined_height = max(height_single * 1.9, 10))
}

for (id in model_ids) {
  cache <- readRDS(file.path(cache_dir, paste0(id, ".rds")))
  message("\n--- ", id, " from cached epred draws ---")

  plots <- switch(
    cache$kind,
    intercept = plot_intercept(cache),
    continuous = plot_continuous(cache),
    categorical = plot_categorical(cache),
    stop("Unsupported cache kind for orchard plot: ", cache$kind)
  )

  save_plot(plots$location, paste0(id, "_orchard_location"),
            width = plots$width, height = plots$height)
  save_plot(plots$scale, paste0(id, "_orchard_scale"),
            width = plots$width, height = plots$height)
  save_plot(plots$combined, paste0(id, "_orchard_combined"),
            width = plots$width, height = plots$combined_height)
}
