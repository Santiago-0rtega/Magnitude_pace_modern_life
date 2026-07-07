# Wrapper around orchaRd:::.safe_lnM_indep() for a single contrast
get_lnM_safe <- function(m1, m2, s1, s2, n1, n2,
                         min_kept = 100000,
                         chunk_init = 5000,
                         seed = 123) {
  out <- orchaRd:::.safe_lnM_indep(
    x1bar = m1, x2bar = m2,
    sd1 = s1, sd2 = s2,
    n1 = n1, n2 = n2,
    min_kept = min_kept,
    chunk_init = chunk_init,
    seed = seed
  )

  tibble::tibble(
    yi_lnM_safe = out$lnM_SAFE,
    vi_lnM_safe = out$var_lnM_SAFE,
    draws_kept  = out$kept,
    draws_total = out$total,
    attempts    = out$attempts,
    status      = out$status
  )
}

# Compute SAFE lnM for every row; cache the result
compute_or_read_lnm_safe <- function(dat_lnm_input,
                                     out_path   = here::here("Rdata", "effect_sizes", "proceed_lnm_safe.rds"),
                                     diag_path  = here::here("Rdata", "tables", "safe_lnm_diagnostics.csv"),
                                     recompute  = FALSE,
                                     min_kept   = 100000,
                                     chunk_init = 5000,
                                     seed       = 123) {
  if (file.exists(out_path) && !recompute) {
    message("Reading cached SAFE lnM: ", out_path)
    return(readRDS(out_path))
  }

  # Read-only book render: never recompute effect sizes on the fly.
  if (isTRUE(getOption("pace.read_only", FALSE)) && !recompute) {
    stop("[read-only render] SAFE lnM cache not found: ", out_path,
         "\nThe book will not recompute effect sizes. Sync the .rds, or set recompute = TRUE deliberately.",
         call. = FALSE)
  }

  message("Computing SAFE lnM for ", nrow(dat_lnm_input), " contrasts ...")

  dat_es <- dat_lnm_input |>
    dplyr::mutate(
      es_id_db = factor(es_id),
      safe = purrr::pmap(
        list(mean1, mean2, sd1, sd2, n1, n2),
        ~ get_lnM_safe(..1, ..2, ..3, ..4, ..5, ..6,
                       min_kept   = min_kept,
                       chunk_init = chunk_init,
                       seed       = seed)
      )
    ) |>
    tidyr::unnest(safe) |>
    dplyr::filter(
      is.finite(yi_lnM_safe),
      is.finite(vi_lnM_safe),
      vi_lnM_safe > 0
    ) |>
    dplyr::mutate(
      es_id_model = factor(seq_len(dplyr::n()))
    )

  # Save diagnostics
  safe_diag <- dat_es |>
    dplyr::select(es_id_db, es_id_model, status, draws_kept, draws_total, attempts)

  readr::write_csv(safe_diag, diag_path)

  saveRDS(dat_es, out_path)
  dat_es
}
