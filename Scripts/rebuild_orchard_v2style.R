# Redraws the time-adjusted (additive) orchard figures in the style of the primary
# orchard figures; these are the figures shown in the book appendix.
# Run from the repository root (paths resolved with here::here()).
suppressMessages({
  library(here)
  library(tidyverse)
  library(brms)
  library(tidybayes)
  library(ggbeeswarm)
  library(patchwork)
})

out_fig_dir <- here::here("outputs", "figures", "publication", "orchard_additive")
dat_es <- readRDS(here::here("outputs", "effect_sizes", "proceed_lnm_safe.rds"))

# ---- Exact visual style from chapters/Rscripts/17_plot_orchard_from_epred_cache.R ----
COL_LOCATION       <- "#0072B2"
COL_LOCATION_LIGHT <- "#88CCEE"
COL_SCALE          <- "#D55E00"
COL_SCALE_LIGHT    <- "#E69F00"
PRECISION_LABEL <- "Effect-size precision (1/SE)"
PRECISION_GUIDE <- ggplot2::guide_legend(
  override.aes = list(shape = 21, fill = "white", colour = "grey40", alpha = 1)
)
CATEGORY_COLS <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
  "#56B4E9", "#6A3D9A", "#999999", "#000000"
)

D_EQ_BREAKS <- c(0.2, 0.5, 0.8, sqrt(2))
LNM_REF <- log(D_EQ_BREAKS[1:3] / sqrt(2))
D_EQ_AXIS <- ggplot2::sec_axis(
  ~ sqrt(2) * exp(.),
  breaks = D_EQ_BREAKS,
  labels = c("0.2", "0.5", "0.8", "1.41"),
  name = expression(Approximate~italic(d)[plain(eq)])
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
      plot.title         = ggplot2::element_text(size = 13, face = "bold")
    )
}

apply_labels <- function(x, lbl_map) {
  if (is.null(lbl_map)) return(x)
  lbl_map <- lbl_map[!is.na(names(lbl_map))]
  ifelse(x %in% names(lbl_map), lbl_map[x], x)
}

# Holds `timevar` at its fitted-data mean so the additive time term does not
# distort the level comparison (no precomputed emmeans cache exists for these
# post-hoc models, so estimates are drawn directly via tidybayes).
get_epred_estimates <- function(fit, moderator, timevar, levels_vec, dpar = NULL) {
  nd <- data.frame(x = levels_vec, stringsAsFactors = FALSE)
  names(nd)[1] <- moderator
  nd[[timevar]]     <- mean(fit$data[[timevar]], na.rm = TRUE)
  nd[["es_id_model"]] <- NA
  nd[["ref_id"]]      <- NA
  nd[["sp_ncbi"]]     <- NA
  nd[["sp_ncbi_canonical"]] <- NA

  draws <- tidybayes::epred_draws(fit, newdata = nd, re_formula = NA,
                                  dpar = !is.null(dpar), ndraws = 1000)
  value_col <- if (is.null(dpar)) ".epred" else dpar

  draws |>
    dplyr::group_by(.data[[moderator]]) |>
    dplyr::summarise(
      estimate = mean(.data[[value_col]]),
      lowerCL  = quantile(.data[[value_col]], 0.025),
      upperCL  = quantile(.data[[value_col]], 0.975),
      .groups  = "drop"
    ) |>
    dplyr::rename(level = dplyr::all_of(moderator))
}

save_plot <- function(p, stem, width, height) {
  ggplot2::ggsave(file.path(out_fig_dir, paste0(stem, ".png")), p,
                  width = width, height = height, dpi = 300, type = "cairo")
  ggplot2::ggsave(file.path(out_fig_dir, paste0(stem, ".pdf")), p,
                  width = width, height = height, device = grDevices::cairo_pdf)
}

build_orchard <- function(fit, model_id, moderator, moderator_label,
                          timevar, label_map = NULL, width = 9) {
  levels_vec <- levels(factor(fit$data[[moderator]]))

  ests     <- get_epred_estimates(fit, moderator, timevar, levels_vec, dpar = NULL)
  sig_ests <- get_epred_estimates(fit, moderator, timevar, levels_vec, dpar = "sigma")

  raw <- dat_es[!is.na(dat_es[[moderator]]) & !is.na(dat_es[[timevar]]), ] |>
    dplyr::mutate(level = as.character(.data[[moderator]]),
                  precision = 1 / sqrt(vi_lnM_safe))

  ests$level     <- apply_labels(ests$level, label_map)
  sig_ests$level <- apply_labels(sig_ests$level, label_map)
  raw$level      <- apply_labels(raw$level, label_map)

  lev_order <- ests$level[order(ests$estimate)]
  ests$level     <- factor(ests$level, levels = lev_order)
  sig_ests$level <- factor(sig_ests$level, levels = lev_order)
  raw$level      <- factor(raw$level, levels = lev_order)

  n_levels <- nlevels(ests$level)
  level_cols <- CATEGORY_COLS[seq_len(n_levels)]
  names(level_cols) <- levels(ests$level)
  y_breaks <- seq_len(n_levels)

  ests$y_index     <- as.numeric(ests$level)
  sig_ests$y_index <- as.numeric(sig_ests$level)
  raw$y_index      <- as.numeric(raw$level)
  cap_half_height <- 0.12

  p_loc <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = LNM_REF, linetype = "dotted",
                        colour = "grey72", linewidth = 0.45) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggbeeswarm::geom_quasirandom(
      data = raw,
      ggplot2::aes(x = yi_lnM_safe, y = y_index,
                   size = precision, colour = level, fill = level),
      alpha = 0.30, shape = 21, groupOnX = FALSE
    ) +
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
      linewidth = 1.25, colour = "grey10"
    ) +
    ggplot2::geom_point(
      data = ests,
      ggplot2::aes(x = estimate, y = y_index),
      size = 3.4, shape = 23, fill = "white", colour = "grey10", stroke = 1.1
    ) +
    ggplot2::scale_colour_manual(values = level_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = level_cols, guide = "none") +
    ggplot2::scale_size_continuous(name = PRECISION_LABEL, range = c(0.4, 4),
                                   guide = PRECISION_GUIDE) +
    ggplot2::scale_x_continuous(sec.axis = D_EQ_AXIS) +
    ggplot2::scale_y_continuous(breaks = y_breaks, labels = lev_order) +
    ggplot2::labs(x = "lnM", y = NULL, title = "A)") +
    theme_orchard()

  raw_sig <- raw |>
    dplyr::left_join(dplyr::select(ests, level, loc_est = estimate), by = "level") |>
    dplyr::mutate(abs_resid = abs(yi_lnM_safe - loc_est))

  p_scl <- ggplot2::ggplot() +
    ggbeeswarm::geom_quasirandom(
      data = raw_sig,
      ggplot2::aes(x = abs_resid, y = y_index,
                   size = precision, colour = level, fill = level),
      alpha = 0.30, shape = 21, groupOnX = FALSE
    ) +
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
    ggplot2::scale_size_continuous(name = PRECISION_LABEL, range = c(0.4, 4),
                                   guide = PRECISION_GUIDE) +
    ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_y_continuous(breaks = y_breaks, labels = lev_order) +
    ggplot2::labs(x = "residual lnM (SD)", y = NULL, title = "B)") +
    theme_orchard()

  height_single <- max(2.5 + n_levels * 0.55, 5)
  p_comb <- patchwork::wrap_plots(
    p_loc + ggplot2::theme(legend.position = "none"), p_scl,
    ncol = 1, guides = "collect"
  ) & ggplot2::theme(legend.position = "bottom")

  save_plot(p_loc,  paste0(model_id, "_orchard_location"), width, height_single)
  save_plot(p_scl,  paste0(model_id, "_orchard_scale"),    width, height_single)
  save_plot(p_comb, paste0(model_id, "_orchard_combined"), width, max(height_single * 1.9, 10))
  invisible(p_comb)
}

label_maps <- list(
  disturbance = c("Climate change" = "Climate change", "Hunt_harv" = "Hunting / harvesting",
                  "Introduction" = "Introduction", "Landscapechange" = "Landscape change",
                  "Other" = "Other (in situ natural variation)", "Pollution" = "Pollution",
                  "Responsetointroductions" = "Response to introductions"),
  design      = c("Allochronic" = "Allochronic", "Synchronic" = "Synchronic"),
  trait_type  = c("behaviour" = "Behaviour", "growth" = "Growth", "otherLH" = "Other life history",
                  "othermorphology" = "Other morphology", "phenology" = "Phenology",
                  "physio" = "Physiology", "response" = "Response (performance ratio)",
                  "size" = "Body size"),
  genphen     = c("Genetic" = "Genetic", "Phenotypic" = "Phenotypic")
)

specs <- tibble::tribble(
  ~model_id,                            ~moderator,    ~moderator_label,        ~timevar,            ~width,
  "disturbance_plus_log10_years",       "disturbance", "Disturbance context",   "log10_years",        13,
  "disturbance_plus_log10_generations", "disturbance", "Disturbance context",   "log10_generations",  13,
  "design_plus_log10_years",            "design",      "Comparison design",    "log10_years",        9,
  "design_plus_log10_generations",      "design",      "Comparison design",    "log10_generations",   9,
  "trait_type_plus_log10_years",        "trait_type",  "Trait type",           "log10_years",        13,
  "trait_type_plus_log10_generations",  "trait_type",  "Trait type",           "log10_generations",  13,
  "genphen_plus_log10_years",           "genphen",     "Phenotypic vs. genetic", "log10_years",        9,
  "genphen_plus_log10_generations",     "genphen",     "Phenotypic vs. genetic", "log10_generations",  9
)

for (i in seq_len(nrow(specs))) {
  row <- specs[i, ]
  cat("== rebuilding (v2 style) ", row$model_id, " ==\n")
  fit <- readRDS(here::here("outputs", "models", paste0(row$model_id, "_ls_additive.rds")))
  build_orchard(fit, row$model_id, row$moderator, row$moderator_label,
               row$timevar, label_map = label_maps[[row$moderator]], width = row$width)
  rm(fit); gc()
}

cat("\nDONE v2 style rebuild.\n")
