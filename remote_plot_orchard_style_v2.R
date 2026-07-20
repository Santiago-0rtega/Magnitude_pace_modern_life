setwd("/home/ortegara/Documents/PACE")
library(here)
suppressPackageStartupMessages({
  source(here::here("R", "00_packages.R"))
  source(here::here("R", "06_model_registry.R"))
  source(here::here("R", "09_model_summaries.R"))
})
library(ggbeeswarm)

out_dir <- here::here("outputs", "figures", "orchard")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat_es <- readRDS(here::here("outputs", "effect_sizes", "proceed_lnm_safe.rds"))

cb_cols <- c("#88CCEE","#CC6677","#DDCC77","#117733","#332288",
             "#AA4499","#44AA99","#999933","#882255","#661100",
             "#6699CC","#888888","#E69F00","#56B4E9","#009E73",
             "#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

# Approximate large-sample conversion: |d| = sqrt(2) * exp(lnM).
d_ref <- c(0.2, 0.5, 0.8)
lnm_ref <- log(d_ref / sqrt(2))
d_axis <- ggplot2::sec_axis(~ sqrt(2) * exp(.),
  breaks = c(d_ref, sqrt(2)), labels = c("0.2", "0.5", "0.8", "1.41"),
  name = "Approximate |d|")
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
  lbl_map <- lbl_map[!is.na(names(lbl_map))]
  ifelse(x %in% names(lbl_map), lbl_map[x], x)
}

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
  out <- tibble::tibble(level = ref_level, estimate = intercept,
                        lowerCL = intercept_lo, upperCL = intercept_hi)
  for (r in rows) {
    lvl <- sub(pat, "", names(beta)[r])
    out <- dplyr::bind_rows(out, tibble::tibble(
      level    = lvl,
      estimate = intercept + beta[r],
      lowerCL  = intercept + lower[r],
      upperCL  = intercept + upper[r]))
  }
  out
}

get_sig_estimates <- function(fit, moderator, levels_vec = NULL) {
  if (is.null(levels_vec)) {
    levels_vec <- levels(factor(fit$data[[moderator]]))
  }
  nd <- data.frame(x = levels_vec, stringsAsFactors = FALSE)
  names(nd)[1] <- moderator
  nd[["es_id_model"]] <- NA; nd[["ref_id"]] <- NA; nd[["sp_ncbi"]] <- NA
  tidybayes::epred_draws(
    fit, newdata = nd, re_formula = NA, dpar = TRUE, ndraws = 1000
  ) |>
    dplyr::group_by(.data[[moderator]]) |>
    dplyr::summarise(
      estimate = median(sigma),
      lowerCL = quantile(sigma, 0.025),
      upperCL = quantile(sigma, 0.975),
      .groups = "drop"
    ) |>
    dplyr::rename(level = dplyr::all_of(moderator))
}

get_pred_interval <- function(fit, moderator, levels_vec) {
  nd <- data.frame(x = levels_vec, stringsAsFactors = FALSE)
  names(nd)[1] <- moderator
  nd[["es_id_model"]] <- NA; nd[["ref_id"]] <- NA; nd[["sp_ncbi"]] <- NA
  pp <- tryCatch(
    tidybayes::add_predicted_draws(nd, fit, re_formula = NA, ndraws = 1000),
    error = function(e) { message("PI failed: ", conditionMessage(e)); NULL })
  if (is.null(pp)) return(NULL)
  pp |>
    dplyr::group_by(.data[[moderator]]) |>
    dplyr::summarise(lowerPR = quantile(.prediction, 0.025),
                     upperPR = quantile(.prediction, 0.975), .groups = "drop") |>
    dplyr::rename(level = dplyr::all_of(moderator))
}

save_plot <- function(p, stem, width = 9, height = 6) {
  cairo_pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width, height = height)
  print(p); dev.off()
  png(file.path(out_dir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 300, type = "cairo")
  print(p); dev.off()
  message("Saved: ", stem)
}

# ---- Model registry with PDF-derived labels -----------------------------------
all_models <- list(
  list(id="m01", file="m01_ls_disturbance", moderator="disturbance",
       label="Disturbance context", type="categorical",
       lbl=c(
         "Climate change"  = "Climate change",
         "Hunt_harv"       = "Hunting / harvesting",
         "Introduction"    = "Introduction",
         "Landscapechange" = "Landscape change",
         "Other"           = "Other (natural variation)",
         "Pollution"       = "Pollution",
         "Responsetointroductions" = "Response to introductions"
       )),
  list(id="m02", file="m02_ls_design", moderator="design",
       label="Comparison design", type="categorical",
       lbl=c(
         "Allochronic" = "Allochronic (same pop., time series)",
         "Synchronic"  = "Synchronic (diverged pops.)"
       )),
  list(id="m03", file="m03_ls_log10_years", moderator="log10_years",
       label="Elapsed time (log10 years)", type="continuous", lbl=NULL),
  list(id="m04", file="m04_ls_log10_generations", moderator="log10_generations",
       label="Elapsed time (log10 generations)", type="continuous", lbl=NULL),
  list(id="m05", file="m05_ls_trait_type", moderator="trait_type",
       label="Trait type", type="categorical",
       lbl=c(
         "behaviour"      = "Behaviour",
         "growth"         = "Growth",
         "otherLH"        = "Other life history",
         "othermorphology"= "Other morphology",
         "phenology"      = "Phenology",
         "physio"         = "Physiology",
         "response"       = "Response (performance ratio)",
         "size"           = "Body size"
       )),
  list(id="m06", file="m06_ls_taxa", moderator="taxa",
       label="Taxonomic group", type="categorical",
       lbl=c(
         "Amphibian"="Amphibian","Annelid"="Annelid","Arthropod"="Arthropod",
         "Bird"="Bird","Fish"="Fish","Mammal"="Mammal","Mollusc"="Mollusc",
         "Plant"="Plant","Reptile"="Reptile"
       )),
  list(id="m07", file="m07_ls_genphen", moderator="genphen",
       label="Phenotypic vs. genetic", type="categorical",
       lbl=c(
         "Genetic"    = "Genetic (common garden / QG)",
         "Phenotypic" = "Phenotypic (wild-measured)"
       )),
  list(id="m08", file="m08_ls_env_change", moderator="env_change",
       label="Environmental-change context", type="categorical",
       lbl=c(
         "novel"   = "Novel (defined start point)",
         "ongoing" = "Ongoing (measured within)"
       )),
  list(id="m09", file="m09_ls_data_type", moderator="data_type",
       label="Data type", type="categorical",
       lbl=c(
         "linear"      = "Linear (1D)",
         "area (2D)"   = "Area (2D)",
         "cube (3D)"   = "Volume / mass (3D)",
         "count"       = "Count",
         "proportion"  = "Proportion",
         "time"        = "Time",
         "date"        = "Date (interval)",
         "temperature" = "Temperature",
         "rate"        = "Rate",
         "ad_ratio"    = "Dimensionless ratio",
         "index"       = "Index (ordinal)",
         "other"       = "Other"
       )),
  list(id="m10", file="m10_ls_transf_data", moderator="transf_data",
       label="Transformation status", type="categorical",
       lbl=c(
         "raw"        = "Raw (untransformed)",
         "ord"        = "Ordination scores",
         "arcsin"     = "Arcsine",
         "arcsin.sqr" = "Angular (arcsin-sqrt)",
         "resid"      = "Residuals / centred",
         "ln"         = "Natural log",
         "log10"      = "Log10"
       )),
  list(id="m11", file="m11_ls_data_scale", moderator="data_scale",
       label="Measurement scale", type="categorical",
       lbl=c(
         "interval" = "Interval (arbitrary zero)",
         "ratio"    = "Ratio (true zero)"
       )),
  list(id="m06b", file="m06b_ls_taxa_v2", moderator="taxa_v2",
       label="Taxonomic group (collapsed)", type="categorical",
       lbl=c(
         "Bird"    = "Bird", "Fish"   = "Fish", "Insect" = "Insect",
         "Mammal"  = "Mammal", "Plant" = "Plant", "Reptile" = "Reptile"
       )),
  list(id="m10b", file="m10b_ls_transf_data_v2", moderator="transf_data_v2",
       label="Transformation status (collapsed)", type="categorical",
       lbl=c(
         "Untransformed" = "Untransformed (raw)",
         "Transformed"   = "Transformed (ln / log10 / arcsin / resid / ord)"
       ))
)

# Restrict plotting to explicitly requested, converged models. When the
# environment variable is unset, retain the historical all-model behaviour.
orchard_model_ids <- strsplit(Sys.getenv("ORCHARD_MODEL_IDS", ""), ",", fixed = TRUE)[[1]]
orchard_model_ids <- trimws(orchard_model_ids[nzchar(orchard_model_ids)])
if (length(orchard_model_ids) > 0) {
  all_models <- Filter(function(m) m$id %in% orchard_model_ids, all_models)
}

# ---- m00: intercept-only overall distribution --------------------------------
m00_path <- here::here("outputs", "models", "m00_ls_intercept_only.rds")
if (file.exists(m00_path) &&
    (length(orchard_model_ids) == 0 || "m00" %in% orchard_model_ids)) {
  message("\n--- m00 (intercept-only) ---")
  fit00 <- readRDS(m00_path)

  raw00 <- dat_es |>
    dplyr::mutate(precision = 1 / sqrt(vi_lnM_safe))

  fe00      <- brms::fixef(fit00)
  int_est   <- fe00["Intercept",    "Estimate"]
  int_lo    <- fe00["Intercept",    "Q2.5"]
  int_hi    <- fe00["Intercept",    "Q97.5"]
  sig_est   <- exp(fe00["sigma_Intercept", "Estimate"])

  vc <- brms::VarCorr(fit00)
  cap <- sprintf("n = %d effect sizes | study SD = %.2f | phylogeny SD = %.2f",
                 nrow(raw00),
                 vc$ref_id$sd["Intercept","Estimate"],
                 vc$sp_ncbi$sd["Intercept","Estimate"])

  # Location panel
  p_m00_loc <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = lnm_ref, linetype = "dotted", colour = "grey72", linewidth = 0.45) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.65) +
    ggbeeswarm::geom_quasirandom(data = raw00,
      ggplot2::aes(x = yi_lnM_safe, y = 1, size = precision),
      alpha = 0.20, shape = 21, fill = "#88CCEE", colour = "#0072B2",
      groupOnX = FALSE) +
    ggplot2::geom_linerange(
      ggplot2::aes(y = 1, xmin = int_lo, xmax = int_hi),
      linewidth = 1.5, colour = "grey15") +
    ggplot2::geom_point(
      ggplot2::aes(x = int_est, y = 1),
      size = 5, shape = 21, fill = "white", colour = "grey10", stroke = 1.2) +
    ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.3, 4)) +
    ggplot2::scale_x_continuous(sec.axis = d_axis) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::labs(x = "Location effect (lnM)", y = NULL, title = "A)") +
    theme_orchard()

  # Scale panel: within-overall residuals vs sigma intercept
  sig_lo  <- exp(fe00["sigma_Intercept", "Q2.5"])
  sig_hi  <- exp(fe00["sigma_Intercept", "Q97.5"])

  raw00_sig <- raw00 |>
    dplyr::mutate(abs_resid = abs(yi_lnM_safe - int_est))

  p_m00_scl <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggbeeswarm::geom_quasirandom(data = raw00_sig,
      ggplot2::aes(x = abs_resid, y = 1, size = precision),
      alpha = 0.24, shape = 21, fill = "#E69F00", colour = "#D55E00",
      groupOnX = FALSE) +
    ggplot2::geom_linerange(
      ggplot2::aes(y = 1, xmin = sig_lo, xmax = sig_hi),
      linewidth = 3, colour = "white") +
    ggplot2::geom_linerange(
      ggplot2::aes(y = 1, xmin = sig_lo, xmax = sig_hi),
      linewidth = 1.5, colour = "grey15") +
    ggplot2::geom_point(
      ggplot2::aes(x = sig_est, y = 1),
      size = 5, shape = 21, fill = "white", colour = "grey10", stroke = 1.2) +
    ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.3, 4)) +
    ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::labs(x = "Residual lnM (SD)", y = NULL, title = "B)") +
    theme_orchard()

  p_m00_comb <- patchwork::wrap_plots(
    p_m00_loc + ggplot2::theme(legend.position = "none"),
    p_m00_scl, ncol = 1
  )

  save_plot(p_m00_loc,  "m00_orchard_location", width = 9, height = 4)
  save_plot(p_m00_scl,  "m00_orchard_scale",    width = 9, height = 4)
  save_plot(p_m00_comb, "m00_orchard_combined", width = 9, height = 8)
}

# ---- Main loop ---------------------------------------------------------------
for (m in all_models) {
  rds_path <- here::here("outputs", "models", paste0(m$file, ".rds"))
  if (!file.exists(rds_path)) {
    message("\nSkipping ", m$id, " — RDS not found")
    next
  }
  message("\n--- ", m$id, " (", m$moderator, ") ---")
  fit <- readRDS(rds_path)

  raw <- dat_es[!is.na(dat_es[[m$moderator]]), ] |>
    dplyr::mutate(
      level     = as.character(.data[[m$moderator]]),
      precision = 1 / sqrt(vi_lnM_safe)
    )

  if (m$type == "categorical") {
    ests     <- get_level_estimates(fit, m$moderator)
    sig_ests <- get_sig_estimates(fit, m$moderator, unique(raw$level))
    pi_df    <- get_pred_interval(fit, m$moderator, unique(raw$level))

    if (!is.null(pi_df)) {
      ests <- dplyr::left_join(ests, pi_df, by = "level")
    } else {
      ests$lowerPR <- NA_real_; ests$upperPR <- NA_real_
    }

    # Apply human-readable labels
    if (!is.null(m$lbl)) {
      ests$level     <- apply_labels(ests$level,     m$lbl)
      sig_ests$level <- apply_labels(sig_ests$level, m$lbl)
      raw$level      <- apply_labels(raw$level,      m$lbl)
    }

    # Order levels by location estimate
    lev_order      <- ests$level[order(ests$estimate)]
    ests$level     <- factor(ests$level,     levels = lev_order)
    sig_ests$level <- factor(sig_ests$level, levels = lev_order)
    raw$level      <- factor(raw$level,      levels = lev_order)

    n_levels <- nlevels(ests$level)
    colors   <- cb_cols[seq_len(n_levels)]
    names(colors) <- levels(ests$level)

    # ---- Location plot (orchaRd style) ----------------------------------------
    p_loc <- ggplot2::ggplot() +
      ggplot2::geom_vline(xintercept = lnm_ref, linetype = "dotted", colour = "grey72", linewidth = 0.45) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.65) +
      ggbeeswarm::geom_quasirandom(data = raw,
        ggplot2::aes(x = yi_lnM_safe, y = level,
                     size = precision, colour = level, fill = level),
        alpha = 0.30, shape = 21, groupOnX = FALSE) +
      { if (!is.null(pi_df))
          ggplot2::geom_linerange(data = ests,
            ggplot2::aes(y = level, xmin = lowerPR, xmax = upperPR),
            linewidth = 0.4, colour = "grey40") } +
      ggplot2::geom_linerange(data = ests,
        ggplot2::aes(y = level, xmin = lowerCL, xmax = upperCL),
        linewidth = 1.5, colour = "grey15") +
      ggplot2::geom_point(data = ests,
        ggplot2::aes(x = estimate, y = level),
        size = 4, shape = 21, fill = "white", colour = "grey10", stroke = 1.2) +
      ggplot2::scale_colour_manual(values = colors, guide = "none") +
      ggplot2::scale_fill_manual(values = colors, guide = "none") +
      ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.4, 4)) +
      ggplot2::scale_x_continuous(sec.axis = d_axis) +
      ggplot2::labs(x = "Location effect (lnM)", y = NULL, title = "A)") +
      theme_orchard()

    # ---- Scale plot with within-group residuals --------------------------------
    # Absolute residuals and sigma are both on the nonnegative residual-SD scale.
    raw_sig <- raw |>
      dplyr::left_join(
        dplyr::select(ests, level, loc_est = estimate),
        by = "level") |>
      dplyr::mutate(
        abs_resid = abs(yi_lnM_safe - loc_est)
      )

    p_scl <- ggplot2::ggplot() +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      ggbeeswarm::geom_quasirandom(data = raw_sig,
        ggplot2::aes(x = abs_resid, y = level,
                     size = precision, colour = level, fill = level),
        alpha = 0.30, shape = 21, groupOnX = FALSE) +
      ggplot2::geom_linerange(data = sig_ests,
        ggplot2::aes(y = level, xmin = lowerCL, xmax = upperCL),
        linewidth = 3, colour = "white") +
      ggplot2::geom_linerange(data = sig_ests,
        ggplot2::aes(y = level, xmin = lowerCL, xmax = upperCL),
        linewidth = 1.5, colour = "grey15") +
      ggplot2::geom_point(data = sig_ests,
        ggplot2::aes(x = estimate, y = level),
        size = 4, shape = 21, fill = "white", colour = "grey10", stroke = 1.2) +
      ggplot2::scale_colour_manual(values = colors, guide = "none") +
      ggplot2::scale_fill_manual(values = colors, guide = "none") +
      ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.4, 4)) +
      ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::labs(x = "Residual lnM (SD)", y = NULL, title = "B)") +
      theme_orchard()

    p_comb <- patchwork::wrap_plots(
      p_loc + ggplot2::theme(legend.position = "none"),
      p_scl, ncol = 1
    )

    h <- 2.5 + n_levels * 0.55
    save_plot(p_loc,  paste0(m$id, "_orchard_location"), width = 9, height = max(h, 5))
    save_plot(p_scl,  paste0(m$id, "_orchard_scale"),    width = 9, height = max(h, 5))
    save_plot(p_comb, paste0(m$id, "_orchard_combined"), width = 9, height = max(h * 1.9, 10))

  } else {
    # ---- Continuous models (m03/m04) ------------------------------------------
    rng <- range(raw[[m$moderator]], na.rm = TRUE)
    nd  <- data.frame(x = seq(rng[1], rng[2], length.out = 100))
    names(nd) <- m$moderator
    nd[["es_id_model"]] <- NA; nd[["ref_id"]] <- NA; nd[["sp_ncbi"]] <- NA

    pred <- tidybayes::add_epred_draws(nd, fit, re_formula = NA, ndraws = 500) |>
      dplyr::group_by(.data[[m$moderator]]) |>
      dplyr::summarise(estimate = median(.epred),
                       lowerCL  = quantile(.epred, 0.025),
                       upperCL  = quantile(.epred, 0.975), .groups = "drop")

    pred_sig <- tidybayes::epred_draws(fit, newdata = nd, dpar = TRUE,
                                       re_formula = NA, ndraws = 500) |>
      dplyr::group_by(.data[[m$moderator]]) |>
      dplyr::summarise(estimate = median(sigma),
                       lowerCL  = quantile(sigma, 0.025),
                       upperCL  = quantile(sigma, 0.975), .groups = "drop")

    # Residuals from the predicted location line for scatter on scale plot
    raw_pred_loc <- approx(pred[[m$moderator]], pred$estimate,
                           xout = raw[[m$moderator]], rule = 2)$y
    raw_cont_sig <- raw |>
      dplyr::mutate(abs_resid = abs(yi_lnM_safe - raw_pred_loc))

    p_loc <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = lnm_ref, linetype = "dotted", colour = "grey72", linewidth = 0.45) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.65) +
      ggplot2::geom_point(data = raw,
        ggplot2::aes(x = .data[[m$moderator]], y = yi_lnM_safe, size = precision),
        alpha = 0.22, shape = 21, fill = "#88CCEE", colour = "#0072B2") +
      ggplot2::geom_ribbon(data = pred,
        ggplot2::aes(x = .data[[m$moderator]], ymin = lowerCL, ymax = upperCL),
        alpha = 0.35, fill = "#0072B2") +
      ggplot2::geom_line(data = pred,
        ggplot2::aes(x = .data[[m$moderator]], y = estimate),
        linewidth = 1.1, colour = "#0072B2") +
      ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.3, 4)) +
      ggplot2::scale_y_continuous(sec.axis = d_axis) +
      ggplot2::labs(x = m$label, y = "lnM", title = "A)") +
      theme_orchard()

    p_scl <- ggplot2::ggplot() +
      ggplot2::geom_point(data = raw_cont_sig,
        ggplot2::aes(x = .data[[m$moderator]], y = abs_resid, size = precision),
        alpha = 0.20, shape = 21, fill = "#88CCEE", colour = "#0072B2") +
      ggplot2::geom_ribbon(data = pred_sig,
        ggplot2::aes(x = .data[[m$moderator]], ymin = lowerCL, ymax = upperCL),
        alpha = 0.35, fill = "grey40") +
      ggplot2::geom_line(data = pred_sig,
        ggplot2::aes(x = .data[[m$moderator]], y = estimate),
        linewidth = 1.1, colour = "grey15") +
      ggplot2::scale_size_continuous(name = "Precision (1/SE)", range = c(0.3, 4)) +
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::labs(x = m$label, y = "Residual lnM (SD)", title = "B)") +
      theme_orchard()

    p_comb <- patchwork::wrap_plots(
      p_loc + ggplot2::theme(legend.position = "none"),
      p_scl, ncol = 1
    )

    save_plot(p_loc,  paste0(m$id, "_orchard_location"), width = 8, height = 5)
    save_plot(p_scl,  paste0(m$id, "_orchard_scale"),    width = 8, height = 4)
    save_plot(p_comb, paste0(m$id, "_orchard_combined"), width = 8, height = 9)
  }
}

message("\nAll orchard plots saved to: ", out_dir)
