# Extract MCMC diagnostics from a brmsfit
extract_diagnostics <- function(fit, model_id, moderator, max_treedepth = 15) {
  draws <- posterior::as_draws_array(fit)
  diag  <- posterior::summarise_draws(
    draws,
    rhat      = posterior::rhat,
    ess_bulk  = posterior::ess_bulk,
    ess_tail  = posterior::ess_tail
  )

  nuts_params <- brms::nuts_params(fit)
  n_divergent <- sum(nuts_params$Value[nuts_params$Parameter == "divergent__"])
  treedepth_vals <- nuts_params$Value[nuts_params$Parameter == "treedepth__"]
  n_treedepth <- sum(treedepth_vals >= max_treedepth, na.rm = TRUE)

  tibble::tibble(
    model_id           = model_id,
    moderator          = moderator,
    max_rhat           = max(diag$rhat, na.rm = TRUE),
    min_bulk_ess       = min(diag$ess_bulk, na.rm = TRUE),
    min_tail_ess       = min(diag$ess_tail, na.rm = TRUE),
    n_divergent        = n_divergent,
    max_treedepth_hits = n_treedepth
  )
}

# Save combined diagnostics table
save_combined_diagnostics <- function(all_diag,
                                      tables_dir = here::here("Rdata", "tables")) {
  tbl <- dplyr::bind_rows(all_diag)
  readr::write_csv(tbl, file.path(tables_dir, "model_diagnostics_all_moderators.csv"))
  tbl
}

# Simple flag: TRUE if diagnostics suggest convergence problems
has_convergence_issues <- function(diag_row,
                                   rhat_thresh = 1.01,
                                   ess_thresh  = 400,
                                   div_thresh  = 10) {
  diag_row$max_rhat > rhat_thresh |
    diag_row$min_bulk_ess < ess_thresh |
    diag_row$min_tail_ess < ess_thresh |
    diag_row$n_divergent > div_thresh
}
