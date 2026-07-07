theme_supplement <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey92"),
      axis.ticks.y = ggplot2::element_blank()
    )
}

# Orchard-like plot for categorical location effects
plot_categorical_location <- function(fe_loc, moderator_label) {
  fe_loc |>
    dplyr::filter(submodel == "location") |>
    ggplot2::ggplot(ggplot2::aes(x = estimate, y = term)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = q2_5, xmax = q97_5),
      size = 0.6
    ) +
    ggplot2::labs(
      x = "Location effect (lnM scale)",
      y = NULL,
      title = paste("Location effects:", moderator_label)
    ) +
    theme_supplement()
}

# Orchard-like plot for categorical scale effects
plot_categorical_scale <- function(fe_scl, moderator_label) {
  fe_scl |>
    dplyr::filter(submodel == "scale") |>
    ggplot2::ggplot(ggplot2::aes(x = estimate, y = term)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = q2_5, xmax = q97_5),
      size = 0.6, colour = "#e06c75"
    ) +
    ggplot2::labs(
      x = "Scale effect (log sigma scale)",
      y = NULL,
      title = paste("Scale effects:", moderator_label),
      caption = "Positive = greater residual heterogeneity. Negative = more consistent divergence."
    ) +
    theme_supplement()
}

# Continuous location prediction ribbon
plot_continuous_location <- function(fit, moderator, moderator_label,
                                     dat_model, n_grid = 100) {
  rng <- range(dat_model[[moderator]], na.rm = TRUE)
  nd  <- tibble::tibble(!!moderator := seq(rng[1], rng[2], length.out = n_grid))
  nd[["es_id_model"]] <- NA
  nd[["ref_id"]]      <- NA
  nd[["sp_ncbi"]]     <- NA

  pred <- tidybayes::add_epred_draws(nd, fit, re_formula = NA, ndraws = 500) |>
    dplyr::group_by(.data[[moderator]]) |>
    dplyr::summarise(
      estimate = median(.epred),
      q2_5     = quantile(.epred, 0.025),
      q97_5    = quantile(.epred, 0.975),
      .groups  = "drop"
    )

  ggplot2::ggplot(pred, ggplot2::aes(x = .data[[moderator]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q2_5, ymax = q97_5), alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = estimate), linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::labs(
      x = moderator_label,
      y = "Predicted lnM",
      title = paste("Location:", moderator_label)
    ) +
    theme_supplement()
}

# Continuous scale prediction ribbon
plot_continuous_scale <- function(fit, moderator, moderator_label,
                                  dat_model, n_grid = 100) {
  rng <- range(dat_model[[moderator]], na.rm = TRUE)
  nd  <- tibble::tibble(!!moderator := seq(rng[1], rng[2], length.out = n_grid))
  nd[["es_id_model"]] <- NA
  nd[["ref_id"]]      <- NA
  nd[["sp_ncbi"]]     <- NA

  pred <- tidybayes::add_epred_draws(nd, fit,
                                     dpar = "sigma",
                                     re_formula = NA,
                                     ndraws = 500) |>
    dplyr::group_by(.data[[moderator]]) |>
    dplyr::summarise(
      estimate = median(.epred),
      q2_5     = quantile(.epred, 0.025),
      q97_5    = quantile(.epred, 0.975),
      .groups  = "drop"
    )

  ggplot2::ggplot(pred, ggplot2::aes(x = .data[[moderator]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q2_5, ymax = q97_5),
                         alpha = 0.25, fill = "#e06c75") +
    ggplot2::geom_line(ggplot2::aes(y = estimate),
                       linewidth = 1, colour = "#e06c75") +
    ggplot2::labs(
      x = moderator_label,
      y = "Predicted sigma (residual SD)",
      title = paste("Scale:", moderator_label),
      caption = "Higher sigma = greater residual heterogeneity."
    ) +
    theme_supplement()
}

# Combined location + scale panel (patchwork)
plot_combined_ls <- function(p_loc, p_scl) {
  patchwork::wrap_plots(p_loc, p_scl, ncol = 1)
}
