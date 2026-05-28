read_proceed <- function(path = here::here("data", "PROCEEDv6.2_RatesDB.csv")) {
  readr::read_csv(path, show_col_types = FALSE) |>
    janitor::clean_names()
}

clean_proceed <- function(dat_raw) {
  dat_raw |>
    dplyr::mutate(
      ref_id      = factor(ref_id),
      sys_id      = factor(sys_id),
      sp_ncbi     = factor(sp_ncbi),
      design      = factor(design),
      genphen     = factor(genphen),
      disturbance = factor(disturbance),
      trait_type  = factor(trait_type),
      taxa        = factor(taxa),
      env_change  = factor(env_change),
      data_type   = factor(data_type),
      transf_data = factor(transf_data),
      data_scale  = factor(data_scale)
    )
}
