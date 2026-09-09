setwd(here::here())
source(here::here("Scripts", "01_paths.R"))
suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
source(here::here("Scripts", "06_model_registry.R"))
source(here::here("Scripts", "10_model_diagnostics.R"))

model_dir <- dir_out("models")
table_dir <- dir_out("tables", "primary")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

model_files <- c(
  disturbance = "m01_ls_disturbance.rds",
  design = "m02_ls_design.rds",
  log10_years = "m03_ls_log10_years.rds",
  log10_generations = "m04_ls_log10_generations.rds",
  trait_type = "m05_ls_trait_type.rds",
  genphen = "m07_ls_genphen.rds",
  env_change = "m08_ls_env_change.rds"
)

fixed_effects <- list()
diagnostics <- list()
for (moderator in names(model_files)) {
  path <- file.path(model_dir, model_files[[moderator]])
  stopifnot(file.exists(path))
  fit <- readRDS(path)
  model_id <- paste0("primary_", moderator)

  fe <- as.data.frame(brms::fixef(fit)) |>
    tibble::rownames_to_column("term") |>
    dplyr::transmute(
      model_id = model_id,
      moderator = moderator,
      submodel = ifelse(grepl("^sigma_", term), "scale", "location"),
      term,
      estimate = Estimate,
      est_error = Est.Error,
      q2_5 = Q2.5,
      q97_5 = Q97.5
    )
  fixed_effects[[moderator]] <- fe
  diagnostics[[moderator]] <- extract_diagnostics(fit, model_id, moderator)
  rm(fit)
  gc(verbose = FALSE)
}

readr::write_csv(dplyr::bind_rows(fixed_effects),
                 file.path(table_dir, "primary_fixef.csv"))
readr::write_csv(dplyr::bind_rows(diagnostics),
                 file.path(table_dir, "primary_diagnostics.csv"))
message("Refreshed matched-primary summaries from ", length(model_files), " current reruns.")
