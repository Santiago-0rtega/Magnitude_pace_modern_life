# Apply lnM inclusion criteria: requires mean1/2, sd1/2, n1/2
require_lnm_vars <- function(dat) {
  dat |>
    dplyr::filter(
      !is.na(mean1), !is.na(mean2),
      !is.na(sd1),   !is.na(sd2),
      !is.na(n1),    !is.na(n2)
    )
}

# Apply the primary generations filter. Eligibility cannot be established when
# elapsed generations are missing, so those rows are excluded.
filter_generations <- function(dat, max_gen = 300) {
  dat |>
    dplyr::filter(!is.na(generations), generations > 0, generations <= max_gen)
}

# Add derived variables needed for modelling
add_derived_vars <- function(dat) {
  dat |>
    dplyr::mutate(
      n_total            = n1 + n2,
      log10_years        = dplyr::if_else(years > 0, log10(years), NA_real_),
      log10_generations  = dplyr::if_else(
        !is.na(generations) & generations > 0,
        log10(generations),
        NA_real_
      )
    )
}

# Filter rows that have valid years (needed for log10_years moderator)
require_years <- function(dat) {
  dat |> dplyr::filter(!is.na(years), years > 0)
}
