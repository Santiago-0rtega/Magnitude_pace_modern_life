# Extract fixed-effects summary from a brmsfit, labelled by model_id and moderator
extract_fixed_effects <- function(fit, model_id, moderator) {
  fe <- brms::fixef(fit, summary = TRUE) |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    tibble::as_tibble() |>
    dplyr::rename(
      estimate   = Estimate,
      est_error  = Est.Error,
      q2_5       = Q2.5,
      q97_5      = Q97.5
    ) |>
    dplyr::mutate(
      model_id  = model_id,
      moderator = moderator,
      submodel  = dplyr::if_else(
        grepl("^sigma_", term),
        "scale",
        "location"
      )
    ) |>
    dplyr::select(model_id, moderator, submodel, term, estimate, est_error, q2_5, q97_5)
  fe
}

# Convenience wrapper: split into location and scale tables
split_location_scale <- function(fe_tbl) {
  list(
    location = dplyr::filter(fe_tbl, submodel == "location"),
    scale    = dplyr::filter(fe_tbl, submodel == "scale")
  )
}

# Accumulate fixed effects across all models and save combined tables
save_combined_summaries <- function(all_fe,
                                    tables_dir = here::here("Rdata", "tables")) {
  loc <- dplyr::bind_rows(lapply(all_fe, `[[`, "location"))
  scl <- dplyr::bind_rows(lapply(all_fe, `[[`, "scale"))

  readr::write_csv(loc, file.path(tables_dir, "location_effects_all_moderators.csv"))
  readr::write_csv(scl, file.path(tables_dir, "scale_effects_all_moderators.csv"))

  list(location = loc, scale = scl)
}

# Extract variance-component summary
extract_variance_components <- function(fit, model_id, moderator) {
  vc <- brms::VarCorr(fit, summary = TRUE)

  purrr::map_dfr(names(vc), function(group) {
    tryCatch({
      sd_arr <- vc[[group]]$sd
      # sd_arr is a named matrix: rows = statistics, cols = random-effect terms
      col_nm <- colnames(sd_arr)[1]
      tibble::tibble(
        model_id  = model_id,
        moderator = moderator,
        group     = group,
        sd_mean   = sd_arr["Estimate", col_nm],
        sd_q2_5   = sd_arr["Q2.5",     col_nm],
        sd_q97_5  = sd_arr["Q97.5",    col_nm]
      )
    }, error = function(e) {
      tibble::tibble(
        model_id  = model_id,
        moderator = moderator,
        group     = group,
        sd_mean   = NA_real_,
        sd_q2_5   = NA_real_,
        sd_q97_5  = NA_real_
      )
    })
  })
}
