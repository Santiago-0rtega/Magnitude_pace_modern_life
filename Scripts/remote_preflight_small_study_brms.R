required <- c("brms", "cmdstanr", "dplyr", "here", "posterior", "readr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
parse(file = "remote_fit_one_small_study_brms.R")
cat("PREFLIGHT_OK\n")
