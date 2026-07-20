parse(file = "remote_rebuild_safe.R")
parse(file = "remote_fit_one_primary.R")
parse(file = "fit_one_sensitivity.R")
source("R/03_filters.R")
x <- data.frame(generations = c(NA, -1, 0, 1, 300, 301))
stopifnot(identical(filter_generations(x)$generations, c(1, 300)))
source("R/06_model_registry.R")
stopifnot(identical(
  moderator_grid$model_id,
  c("m01", "m02", "m03", "m04", "m05", "m07", "m08", "m11")
))
cat("VALIDATION_OK\n")
