# Small-study-effects diagnostic for SAFE lnM.
#
# Following Nakagawa et al. (2022) and the lnM worked example, the predictor is
# based only on group sample sizes. This avoids using an SE that shares
# estimated variance components with lnM itself.

prepare_small_study_data <- function(
    effect_size_path = here::here("Rdata", "effect_sizes", "proceed_lnm_safe.rds"),
    A_path = here::here("Rdata", "phylogeny", "proceed_A_matrix.rds"),
    name_map_path = here::here("Rdata", "phylogeny", "proceed_name_map.rds"),
    phylogeny_script = NULL) {

  if (is.null(phylogeny_script)) {
    candidates <- c(
      here::here("chapters", "Rscripts", "05_phylogeny.R"),
      here::here("R", "05_phylogeny.R")
    )
    phylogeny_script <- candidates[file.exists(candidates)][1]
  }
  if (is.na(phylogeny_script) || !file.exists(phylogeny_script)) {
    stop("Cannot find 05_phylogeny.R in the local or remote project layout.")
  }
  source(phylogeny_script)

  dat <- readRDS(effect_size_path)
  A_full <- readRDS(A_path)
  name_map <- readRDS(name_map_path)

  dat <- apply_phylo_name_map(dat, name_map) |>
    dplyr::filter(
      is.finite(yi_lnM_safe),
      is.finite(vi_lnM_safe), vi_lnM_safe > 0,
      is.finite(n1), n1 > 0,
      is.finite(n2), n2 > 0,
      !is.na(sp_ncbi_canonical), nzchar(sp_ncbi_canonical)
    ) |>
    dplyr::mutate(
      n0 = (n1 * n2) / (n1 + n2),
      n_se = 1 / sqrt(n0),
      n_v = 1 / n0,
      ref_id = droplevels(factor(ref_id)),
      sp_ncbi = factor(sp_ncbi_canonical),
      es_id_model = factor(seq_len(dplyr::n()))
    )

  species <- levels(dat$sp_ncbi)
  A <- expand_A_with_unmatched(A_full, species)

  list(data = dat, A = A)
}

fit_small_study_model <- function(dat, A, predictor = c("n_se", "n_v")) {
  predictor <- match.arg(predictor)

  metafor::rma.mv(
    yi = yi_lnM_safe,
    V = vi_lnM_safe,
    mods = stats::as.formula(paste("~", predictor)),
    random = list(
      ~ 1 | ref_id,
      ~ 1 | sp_ncbi,
      ~ 1 | es_id_model
    ),
    R = list(sp_ncbi = A),
    Rscale = "cor",
    method = "REML",
    test = "t",
    sparse = TRUE,
    data = dat
  )
}

tidy_small_study_model <- function(model, model_name) {
  coefs <- data.frame(
    term = rownames(model$beta),
    estimate = as.numeric(model$beta),
    std_error = as.numeric(model$se),
    ci_lower = as.numeric(model$ci.lb),
    ci_upper = as.numeric(model$ci.ub),
    p_value = as.numeric(model$pval),
    row.names = NULL
  )

  dplyr::mutate(coefs, model = model_name, .before = 1)
}

build_small_study_figure <- function(model, dat,
                                     predictor = c("n_se", "n_v"),
                                     bubble_fill = "#DC143C",
                                     line_colour = bubble_fill) {
  predictor <- match.arg(predictor)
  grid <- data.frame(x = seq(0, max(dat[[predictor]]), length.out = 200))
  pred <- stats::predict(model, newmods = grid$x)
  grid$estimate <- as.numeric(pred$pred)
  grid$ci_lower <- as.numeric(pred$ci.lb)
  grid$ci_upper <- as.numeric(pred$ci.ub)

  x_label <- if (predictor == "n_se") {
    expression(1 / sqrt(n[0]))
  } else {
    expression(1 / n[0])
  }

  ggplot2::ggplot(dat, ggplot2::aes(x = .data[[predictor]], y = yi_lnM_safe)) +
    ggplot2::geom_point(
      ggplot2::aes(size = 1 / sqrt(vi_lnM_safe)),
      shape = 21, fill = bubble_fill, colour = line_colour,
      alpha = 0.16, stroke = 0.2
    ) +
    ggplot2::geom_ribbon(
      data = grid,
      ggplot2::aes(x = x, y = estimate, ymin = ci_lower, ymax = ci_upper),
      inherit.aes = FALSE, fill = bubble_fill, alpha = 0.18
    ) +
    ggplot2::geom_line(
      data = grid, ggplot2::aes(x = x, y = estimate),
      inherit.aes = FALSE, colour = line_colour, linewidth = 1
    ) +
    ggplot2::scale_size_continuous(
      name = "Precision (1/SE)", range = c(0.5, 3.5),
      trans = "sqrt",
      guide = ggplot2::guide_legend(
        override.aes = list(alpha = 1, fill = NA, colour = "black")
      )
    ) +
    ggplot2::labs(
      x = x_label,
      y = "lnM",
      title = NULL, subtitle = NULL, caption = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(legend.position = "right")
}

run_small_study_analysis <- function(
    table_path = here::here("Rdata", "tables", "small_study_effects.csv"),
    model_path = here::here("Rdata", "models", "small_study_effects.rds"),
    figure_path = here::here("Rdata", "figures", "publication", "small_study_effects.png")) {

  inputs <- prepare_small_study_data()
  fit_se <- fit_small_study_model(inputs$data, inputs$A, "n_se")
  fit_v <- fit_small_study_model(inputs$data, inputs$A, "n_v")

  results <- dplyr::bind_rows(
    tidy_small_study_model(fit_se, "Small-study slope: 1 / sqrt(n0)"),
    tidy_small_study_model(fit_v, "Adjustment model: 1 / n0")
  )

  dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(model_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(results, table_path)
  saveRDS(list(n = nrow(inputs$data), fit_se = fit_se, fit_v = fit_v), model_path)

  p <- build_small_study_figure(fit_se, inputs$data, predictor = "n_se")
  ggplot2::ggsave(figure_path, p, width = 8.5, height = 6, dpi = 300)

  invisible(list(data = inputs$data, fit_se = fit_se, fit_v = fit_v,
                 results = results, figure = p))
}

# ══════════════════════════════════════════════════════════════════════════════
# brms location–scale consumers
# ══════════════════════════════════════════════════════════════════════════════
# The primary small-study models are the brms location–scale fits produced on
# totoro by remote_fit_one_small_study_brms.R. Chapters read the small draw
# caches built by remote_precompute_small_study_brms_cache.R and never load the
# ~400 MB fits. The metafor functions above are retained only for the
# frequentist cross-check reported alongside them; rma.mv is homoscedastic, so
# it reproduces the `sigma ~ 1` specification.

SMALL_STUDY_COL <- c(n_se = "#DC143C", n_v = "#FF8C00")

small_study_x_label <- function(moderator) {
  if (moderator == "n_se") expression(1 / sqrt(n[0])) else expression(1 / n[0])
}

read_small_study_cache <- function(id,
                                   cache_dir = here::here("Rdata", "epred_draws")) {
  f <- file.path(cache_dir, paste0("small_study_", id, ".rds"))
  if (!file.exists(f)) return(NULL)
  cache <- readRDS(f)
  if (!identical(cache$summary_spec$version, 2L)) {
    stop("Legacy small-study cache rejected for ", id,
         ". Rebuild with remote_precompute_small_study_brms_cache.R.")
  }
  cache
}

# Posterior fixed effects across one or more caches, with readable term names.
small_study_coef_table <- function(caches) {
  term_labels <- c(
    Intercept       = "Intercept (location)",
    n_se            = "1 / sqrt(n0) slope (location)",
    n_v             = "1 / n0 slope (location)",
    sigma_Intercept = "Intercept (log sigma)",
    sigma_n_se      = "1 / sqrt(n0) slope (log sigma)",
    sigma_n_v       = "1 / n0 slope (log sigma)"
  )

  dplyr::bind_rows(lapply(caches, function(cache) {
    fe <- as.data.frame(cache$fixef)
    nice <- unname(term_labels[rownames(fe)])
    data.frame(
      Model            = cache$id,
      Term             = ifelse(is.na(nice), rownames(fe), nice),
      Estimate         = fe$Estimate,
      `Est. error`     = fe$Est.Error,
      `95% CrI lower`  = fe$Q2.5,
      `95% CrI upper`  = fe$Q97.5,
      row.names        = NULL,
      check.names      = FALSE
    )
  }))
}

# Location panel: raw bubbles sized by precision, with the posterior mean and
# credible bands extended to the predictor's zero bound.
build_small_study_location_panel <- function(cache, show_legend = TRUE) {
  mod <- cache$moderator
  col <- unname(SMALL_STUDY_COL[[mod]])

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_point(
      data = cache$raw,
      ggplot2::aes(x = .data[[mod]], y = yi_lnM_safe,
                   size = 1 / sqrt(vi_lnM_safe)),
      shape = 21, fill = col, colour = col, alpha = 0.16, stroke = 0.2
    ) +
    # show.legend = FALSE: the band levels are stated in the caption. Left on,
    # stat_lineribbon adds a "level" key per panel that guides = "collect"
    # cannot merge across differently coloured panels.
    ggdist::stat_lineribbon(
      data = cache$loc, ggplot2::aes(x = .data[[mod]], y = .epred),
      .width = c(0.50, 0.80, 0.95), alpha = 0.30, fill = col, colour = NA,
      show.legend = FALSE
    ) +
    ggdist::stat_lineribbon(
      data = cache$loc, ggplot2::aes(x = .data[[mod]], y = .epred),
      .width = 0, colour = col, linewidth = 1, show.legend = FALSE
    ) +
    ggplot2::scale_size_continuous(
      name = "Precision (1/SE)", range = c(0.5, 3.5), trans = "sqrt",
      guide = ggplot2::guide_legend(
        override.aes = list(alpha = 1, fill = NA, colour = "black")
      )
    ) +
    ggplot2::labs(x = small_study_x_label(mod), y = "lnM") +
    ggplot2::theme_classic(base_size = 13)

  if (show_legend) p else p + ggplot2::theme(legend.position = "none")
}

# Scale panel: posterior residual SD across the predictor. Flat by construction
# for a `sigma ~ 1` model; informative only for the `sigma ~ x` fits.
build_small_study_scale_panel <- function(cache) {
  mod <- cache$moderator
  col <- unname(SMALL_STUDY_COL[[mod]])

  ggplot2::ggplot(cache$scl, ggplot2::aes(x = .data[[mod]], y = sigma)) +
    ggdist::stat_lineribbon(.width = c(0.50, 0.80, 0.95), alpha = 0.30,
                            fill = col, colour = NA, show.legend = FALSE) +
    ggdist::stat_lineribbon(.width = 0, colour = col, linewidth = 1,
                            show.legend = FALSE) +
    ggplot2::labs(x = small_study_x_label(mod),
                  y = expression("Residual SD (" * sigma * ")")) +
    ggplot2::theme_classic(base_size = 13)
}

# Two-panel publication figure: the primary diagnostic slope beside the
# infinite-sample adjustment model.
build_small_study_brms_figure <- function(cache_se, cache_v) {
  p_se <- build_small_study_location_panel(cache_se, show_legend = TRUE)
  p_v  <- build_small_study_location_panel(cache_v, show_legend = FALSE)

  ((p_se + p_v) / patchwork::guide_area()) +
    patchwork::plot_layout(heights = c(1, 0.10), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "center",
      plot.tag = ggplot2::element_text(face = "bold", size = 14),
      plot.tag.position = c(0.01, 0.99)
    )
}
