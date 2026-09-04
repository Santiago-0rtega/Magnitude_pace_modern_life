invisible(parse(file = "remote_rebuild_safe.R"))
invisible(parse(file = "remote_fit_one_primary.R"))
invisible(parse(file = "remote_fit_one_sensitivity.R"))
source(here::here("chapters", "Rscripts", "03_filters.R"))
x <- data.frame(generations = c(NA, -1, 0, 1, 300, 301))
stopifnot(identical(filter_generations(x)$generations, c(1, 300)))
source(here::here("chapters", "Rscripts", "06_model_registry.R"))

# Primary reported grid. m06/m08/m09/m10 sit in excluded_grid with their
# rationale; m08 is still fitted and rendered as an appendix chapter.
stopifnot(identical(
  moderator_grid$model_id,
  c("m01", "m02", "m03", "m04", "m05", "m07")
))
stopifnot(identical(
  excluded_grid$model_id,
  c("m06", "m08", "m09", "m10")
))
cat("VALIDATION_OK\n")
