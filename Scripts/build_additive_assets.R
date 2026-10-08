# Builds summaries (additive_summaries.rds), diagnostics (additive_diagnostics.csv) and
# initial orchard figures for the time-adjusted (additive) models.
# Run from the repository root (paths resolved with here::here()).
suppressMessages({
  library(here)
  library(tidyverse)
  library(brms)
  library(tidybayes)
  library(ggbeeswarm)
  library(patchwork)
})
source(here::here("Scripts", "10_model_diagnostics.R"))

out_fig_dir   <- here::here("outputs", "figures", "publication", "orchard_additive")
out_tab_dir   <- here::here("outputs", "tables")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)

dat_es <- readRDS(here::here("outputs", "effect_sizes", "proceed_lnm_safe.rds"))

cb_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
             "#AA4499", "#44AA99", "#999933")

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

apply_labels <- function(x, lbl_map) {
  if (is.null(lbl_map)) return(x)
  lbl_map <- lbl_map[!is.na(names(lbl_map))]
  ifelse(x %in% names(lbl_map), lbl_map[x], x)
}

# Unified epred-based estimator: holds `timevar` at its fitted-data mean so the
# additive time term does not distort the level comparison.
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

get_pred_interval <- function(fit, moderator, timevar, levels_vec) {
  nd <- data.frame(x = levels_vec, stringsAsFactors = FALSE)
  names(nd)[1] <- moderator
  nd[[timevar]]     <- mean(fit$data[[timevar]], na.rm = TRUE)
  nd[["es_id_model"]] <- NA
  nd[["ref_id"]]      <- NA
  nd[["sp_ncbi"]]     <- NA
  nd[["sp_ncbi_canonical"]] <- NA

  pp <- tidybayes::add_predicted_draws(nd, fit, re_formula = NA, ndraws = 1000)
  pp |>
    dplyr::group_by(.data[[moderator]]) |>
    dplyr::summarise(
      lowerPR = quantile(.prediction, 0.025),
      upperPR = quantile(.prediction, 0.975),
      .groups = "drop"
    ) |>
    dplyr::rename(level = dplyr::all_of(moderator))
}

build_orchard <- function(fit, model_id, moderator, moderator_label,
                          timevar, label_map = NULL) {
  levels_vec <- levels(factor(fit$data[[moderator]]))

  ests     <- get_epred_estimates(fit, moderator, timevar, levels_vec, dpar = NULL)
  sig_ests <- get_epred_estimates(fit, moderator, timevar, levels_vec, dpar = "sigma")
  pi_df    <- get_pred_interval(fit, moderator, timevar, levels_vec)
  ests     <- dplyr::left_join(ests, pi_df, by = "level")

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
  colors <- cb_cols[seq_len(n_levels)]
  names(colors) <- levels(ests$level)

  p_loc <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggbeeswarm::geom_quasirandom(
      data = raw,
      ggplot2::aes(x = yi_lnM_safe, y = level,
                   size = precision, colour = level, fill = level),
      alpha = 0.30, shape = 21, groupOnX = FALSE
    ) +
    ggplot2::geom_linerange(
      data = ests, ggplot2::aes(y = level, xmin = lowerPR, xmax = upperPR),
      linewidth = 0.4, colour = "grey40"
    ) +
    ggplot2::geom_linerange(
      data = ests, ggplot2::aes(y = level, xmin = lowerCL, xmax = upperCL),
      linewidth = 1.5, colour = "grey15"
    ) +
    ggplot2::geom_point(
      data = ests, ggplot2::aes(x = estimate, y = level),
      size = 4, shape = 21, fill = "white", colour = "grey10", stroke = 1.2
    ) +
    ggplot2::scale_colour_manual(values = colors, guide = "none") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.4, 4)) +
    ggplot2::labs(
      x = "Location effect (lnM), time held at its mean",
      y = NULL, title = paste("Location --", moderator_label)
    ) +
    theme_orchard()

  raw_sig <- raw |>
    dplyr::left_join(dplyr::select(ests, level, loc_est = estimate), by = "level") |>
    dplyr::mutate(abs_resid = abs(yi_lnM_safe - loc_est))

  p_scl <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggbeeswarm::geom_quasirandom(
      data = raw_sig,
      ggplot2::aes(x = abs_resid, y = level,
                   size = precision, colour = level, fill = level),
      alpha = 0.30, shape = 21, groupOnX = FALSE
    ) +
    ggplot2::geom_linerange(
      data = sig_ests, ggplot2::aes(y = level, xmin = lowerCL, xmax = upperCL),
      linewidth = 3, colour = "white"
    ) +
    ggplot2::geom_linerange(
      data = sig_ests, ggplot2::aes(y = level, xmin = lowerCL, xmax = upperCL),
      linewidth = 1.5, colour = "grey15"
    ) +
    ggplot2::geom_point(
      data = sig_ests, ggplot2::aes(x = estimate, y = level),
      size = 4, shape = 21, fill = "white", colour = "grey10", stroke = 1.2
    ) +
    ggplot2::scale_colour_manual(values = colors, guide = "none") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.4, 4)) +
    ggplot2::labs(
      x = "Residual lnM (SD)",
      y = NULL, title = paste("Scale --", moderator_label),
      caption = "Bubbles: absolute residual lnM, sized by precision. Trunk = predicted residual SD (sigma), time held at its mean."
    ) +
    theme_orchard()

  p_comb <- patchwork::wrap_plots(p_loc, p_scl, ncol = 1)
  height_single <- max(2.5 + n_levels * 0.55, 5)

  ggplot2::ggsave(file.path(out_fig_dir, paste0(model_id, "_orchard_combined.png")),
                  p_comb, width = 9, height = max(height_single * 1.9, 10),
                  dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(out_fig_dir, paste0(model_id, "_orchard_combined.pdf")),
                  p_comb, width = 9, height = max(height_single * 1.9, 10),
                  device = cairo_pdf)
  invisible(p_comb)
}

specs <- tibble::tribble(
  ~model_id,                            ~moderator,    ~moderator_label,        ~timevar,
  "disturbance_plus_log10_years",       "disturbance", "Disturbance context",   "log10_years",
  "disturbance_plus_log10_generations", "disturbance", "Disturbance context",   "log10_generations",
  "design_plus_log10_years",            "design",      "Comparison design",    "log10_years",
  "design_plus_log10_generations",      "design",      "Comparison design",    "log10_generations",
  "trait_type_plus_log10_years",        "trait_type",  "Trait type",           "log10_years",
  "trait_type_plus_log10_generations",  "trait_type",  "Trait type",           "log10_generations",
  "genphen_plus_log10_years",           "genphen",     "Phenotypic vs. genetic", "log10_years",
  "genphen_plus_log10_generations",     "genphen",     "Phenotypic vs. genetic", "log10_generations"
)

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

diag_rows <- list()
summary_txts <- list()

for (i in seq_len(nrow(specs))) {
  row <- specs[i, ]
  model_path <- here::here("outputs", "models", paste0(row$model_id, "_ls_additive.rds"))
  cat("== ", row$model_id, " ==\n")
  fit <- readRDS(model_path)

  summary_txts[[row$model_id]] <- capture.output(print(summary(fit)))

  d <- extract_diagnostics(fit, row$model_id, row$moderator)
  d$timevar <- row$timevar
  diag_rows[[row$model_id]] <- d

  build_orchard(fit, row$model_id, row$moderator, row$moderator_label,
               row$timevar, label_map = label_maps[[row$moderator]])

  rm(fit); gc()
}

diag_tbl <- dplyr::bind_rows(diag_rows) |>
  dplyr::select(model_id, moderator, timevar, max_rhat, min_bulk_ess,
               min_tail_ess, n_divergent, max_treedepth_hits)

readr::write_csv(diag_tbl, file.path(out_tab_dir, "additive_diagnostics.csv"))
saveRDS(summary_txts, file.path(out_tab_dir, "additive_summaries.rds"))

cat("\nDONE. Diagnostics:\n")
print(as.data.frame(diag_tbl), row.names = FALSE)
