setwd(here::here())
source(here::here("chapters", "Rscripts", "01_paths.R"))
library(here)
suppressPackageStartupMessages({
  source(here::here("chapters", "Rscripts", "00_packages_fit.R"))
  source(here::here("chapters", "Rscripts", "06_model_registry.R"))
  source(here::here("chapters", "Rscripts", "09_model_summaries.R"))
})
library(ggbeeswarm)

out_dir <- dir_out("figures", "orchard")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat_es <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))

# Colorblind palette (orchaRd default)
cb_cols <- c("#88CCEE","#CC6677","#DDCC77","#117733","#332288",
             "#AA4499","#44AA99","#999933","#882255","#661100",
             "#6699CC","#888888","#E69F00","#56B4E9","#009E73",
             "#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

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

# Marginal means + CI per level from brms fixed effects (treatment coding)
get_level_estimates <- function(fit, moderator) {
  fe    <- brms::fixef(fit)
  beta  <- fe[, "Estimate"]
  lower <- fe[, "Q2.5"]
  upper <- fe[, "Q97.5"]

  intercept    <- beta["Intercept"]
  intercept_lo <- lower["Intercept"]
  intercept_hi <- upper["Intercept"]

  pat  <- paste0("^", moderator)
  rows <- grep(pat, names(beta))

  ref_level <- levels(factor(fit$data[[moderator]]))[1]

  levels_out <- tibble::tibble(
    level    = ref_level,
    estimate = intercept,
    lowerCL  = intercept_lo,
    upperCL  = intercept_hi
  )

  for (r in rows) {
    lvl <- sub(pat, "", names(beta)[r])
    levels_out <- dplyr::bind_rows(levels_out, tibble::tibble(
      level    = lvl,
      estimate = intercept + beta[r],
      lowerCL  = intercept + lower[r],
      upperCL  = intercept + upper[r]
    ))
  }
  levels_out
}

# Prediction interval from posterior_predict (marginalizing RE)
get_pred_interval <- function(fit, moderator, levels_vec) {
  nd <- data.frame(level = levels_vec, stringsAsFactors = FALSE)
  names(nd)[1] <- moderator
  nd[["es_id_model"]] <- NA
  nd[["ref_id"]]      <- NA
  nd[["sp_ncbi"]]     <- NA

  pp <- tryCatch(
    tidybayes::add_predicted_draws(nd, fit, re_formula = NA, ndraws = 1000),
    error = function(e) { message("PI failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(pp)) return(NULL)

  pp |>
    dplyr::group_by(.data[[moderator]]) |>
    dplyr::summarise(
      lowerPR = quantile(.prediction, 0.025),
      upperPR = quantile(.prediction, 0.975),
      .groups = "drop"
    ) |>
    dplyr::rename(level = dplyr::all_of(moderator))
}

save_plot <- function(p, stem, width = 8, height = 6) {
  cairo_pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width, height = height)
  print(p); dev.off()
  png(file.path(out_dir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 300, type = "cairo")
  print(p); dev.off()
  message("Saved: ", stem)
}

# ── Categorical models ─────────────────────────────────────────────────────────
cat_models <- list(
  list(file = "m01_ls_disturbance",  moderator = "disturbance", label = "Disturbance context", id = "m01"),
  list(file = "m02_ls_design",       moderator = "design",      label = "Comparison design",   id = "m02"),
  list(file = "m05_ls_trait_type",   moderator = "trait_type",  label = "Trait type",          id = "m05")
)

for (m in cat_models) {
  message("\n--- ", m$id, " (", m$moderator, ") ---")
  fit <- readRDS(dir_out("models", paste0(m$file, ".rds")))

  raw <- dat_es[!is.na(dat_es[[m$moderator]]), ] |>
    dplyr::mutate(
      level     = as.character(.data[[m$moderator]]),
      precision = 1 / sqrt(vi_lnM_safe)
    )

  ests   <- get_level_estimates(fit, m$moderator)
  pi_df  <- get_pred_interval(fit, m$moderator, unique(raw$level))

  if (!is.null(pi_df)) {
    ests <- dplyr::left_join(ests, pi_df, by = "level")
  } else {
    ests$lowerPR <- NA_real_
    ests$upperPR <- NA_real_
  }

  n_levels <- length(unique(raw$level))
  colors   <- cb_cols[seq_len(n_levels)]
  names(colors) <- sort(unique(raw$level))

  # Location plot — bubbles drawn first (background), then CI lines, then trunk
  p_loc <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggbeeswarm::geom_quasirandom(data = raw,
      ggplot2::aes(x = yi_lnM_safe, y = level,
                   size = precision, colour = level, fill = level),
      alpha = 0.35, shape = 21, groupOnX = FALSE) +
    { if (!is.null(pi_df))
        ggplot2::geom_linerange(data = ests,
          ggplot2::aes(y = level, xmin = lowerPR, xmax = upperPR),
          linewidth = 0.4, colour = "grey40")
    } +
    ggplot2::geom_linerange(data = ests,
      ggplot2::aes(y = level, xmin = lowerCL, xmax = upperCL),
      linewidth = 1.4, colour = "grey20") +
    ggplot2::geom_point(data = ests,
      ggplot2::aes(x = estimate, y = level),
      size = 4, shape = 21, fill = "white", colour = "grey10", stroke = 1.2) +
    ggplot2::scale_colour_manual(values = colors, guide = "none") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.5, 4)) +
    ggplot2::labs(x = "Location effect (lnM)", y = NULL,
                  title = paste("Location —", m$label)) +
    theme_orchard()

  # Scale plot — sigma fixed effects
  fe_all   <- brms::fixef(fit)
  sig_rows <- grep("^sigma_", rownames(fe_all))
  ref_lev  <- levels(factor(fit$data[[m$moderator]]))[1]

  sig_ests <- tibble::tibble(
    level    = ref_lev,
    estimate = fe_all["sigma_Intercept", "Estimate"],
    lowerCL  = fe_all["sigma_Intercept", "Q2.5"],
    upperCL  = fe_all["sigma_Intercept", "Q97.5"]
  )
  for (r in sig_rows) {
    nm <- rownames(fe_all)[r]
    if (nm == "sigma_Intercept") next
    lvl <- sub(paste0("sigma_", m$moderator), "", nm)
    sig_ests <- dplyr::bind_rows(sig_ests, tibble::tibble(
      level    = lvl,
      estimate = fe_all["sigma_Intercept", "Estimate"] + fe_all[nm, "Estimate"],
      lowerCL  = fe_all["sigma_Intercept", "Estimate"] + fe_all[nm, "Q2.5"],
      upperCL  = fe_all["sigma_Intercept", "Estimate"] + fe_all[nm, "Q97.5"]
    ))
  }

  p_scl <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = fe_all["sigma_Intercept", "Estimate"],
                        linetype = "dashed", colour = "grey50") +
    ggplot2::geom_linerange(data = sig_ests,
      ggplot2::aes(y = level, xmin = lowerCL, xmax = upperCL),
      linewidth = 1.4, colour = "#e06c75") +
    ggplot2::geom_point(data = sig_ests,
      ggplot2::aes(x = estimate, y = level),
      size = 4, shape = 21, fill = "white", colour = "#c0392b", stroke = 1.2) +
    ggplot2::labs(x = "Scale effect (log σ)", y = NULL,
                  title = paste("Scale —", m$label),
                  caption = "Positive = greater residual heterogeneity") +
    theme_orchard()

  p_comb <- patchwork::wrap_plots(p_loc, p_scl, ncol = 1)

  h <- 2 + n_levels * 0.55
  save_plot(p_loc,  paste0(m$id, "_orchard_location"), width = 9, height = max(h, 5))
  save_plot(p_scl,  paste0(m$id, "_orchard_scale"),    width = 9, height = max(h * 0.7, 4))
  save_plot(p_comb, paste0(m$id, "_orchard_combined"), width = 9, height = max(h * 1.7, 9))
}

# ── Continuous models ──────────────────────────────────────────────────────────
cont_models <- list(
  list(file = "m03_ls_log10_years",       moderator = "log10_years",
       label = "Elapsed time (log10 years)",       id = "m03"),
  list(file = "m04_ls_log10_generations", moderator = "log10_generations",
       label = "Elapsed time (log10 generations)", id = "m04")
)

for (m in cont_models) {
  message("\n--- ", m$id, " (", m$moderator, ") ---")
  fit <- readRDS(dir_out("models", paste0(m$file, ".rds")))

  raw <- dat_es[!is.na(dat_es[[m$moderator]]), ] |>
    dplyr::mutate(precision = 1 / sqrt(vi_lnM_safe))

  rng <- range(raw[[m$moderator]], na.rm = TRUE)
  nd  <- data.frame(x = seq(rng[1], rng[2], length.out = 100))
  names(nd) <- m$moderator
  nd[["es_id_model"]] <- NA
  nd[["ref_id"]]      <- NA
  nd[["sp_ncbi"]]     <- NA

  pred <- tidybayes::add_epred_draws(nd, fit, re_formula = NA, ndraws = 500) |>
    dplyr::group_by(.data[[m$moderator]]) |>
    dplyr::summarise(
      estimate = median(.epred),
      lowerCL  = quantile(.epred, 0.025),
      upperCL  = quantile(.epred, 0.975),
      .groups  = "drop"
    )

  # Location: bubbles in background, ribbon + line on top
  p_loc <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(data = raw,
      ggplot2::aes(x = .data[[m$moderator]], y = yi_lnM_safe, size = precision),
      alpha = 0.25, shape = 21, fill = "#88CCEE", colour = "#0072B2") +
    ggplot2::geom_ribbon(data = pred,
      ggplot2::aes(x = .data[[m$moderator]], ymin = lowerCL, ymax = upperCL),
      alpha = 0.35, fill = "#0072B2") +
    ggplot2::geom_line(data = pred,
      ggplot2::aes(x = .data[[m$moderator]], y = estimate),
      linewidth = 1.1, colour = "#0072B2") +
    ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.3, 4)) +
    ggplot2::labs(x = m$label, y = "lnM",
                  title = paste("Location —", m$label)) +
    theme_orchard()

  # Scale: sigma ribbon
  pred_sig <- tidybayes::add_epred_draws(nd, fit, dpar = "sigma",
                                         re_formula = NA, ndraws = 500) |>
    dplyr::group_by(.data[[m$moderator]]) |>
    dplyr::summarise(
      estimate = median(.epred),
      lowerCL  = quantile(.epred, 0.025),
      upperCL  = quantile(.epred, 0.975),
      .groups  = "drop"
    )

  p_scl <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = pred_sig,
      ggplot2::aes(x = .data[[m$moderator]], ymin = lowerCL, ymax = upperCL),
      alpha = 0.35, fill = "#e06c75") +
    ggplot2::geom_line(data = pred_sig,
      ggplot2::aes(x = .data[[m$moderator]], y = estimate),
      linewidth = 1.1, colour = "#c0392b") +
    ggplot2::labs(x = m$label, y = "Predicted σ (residual SD)",
                  title = paste("Scale —", m$label),
                  caption = "Higher sigma = greater residual heterogeneity") +
    theme_orchard()

  p_comb <- patchwork::wrap_plots(p_loc, p_scl, ncol = 1)

  save_plot(p_loc,  paste0(m$id, "_orchard_location"), width = 8, height = 5)
  save_plot(p_scl,  paste0(m$id, "_orchard_scale"),    width = 8, height = 4)
  save_plot(p_comb, paste0(m$id, "_orchard_combined"), width = 8, height = 9)
}

message("\nAll orchard plots saved to: ", out_dir)
