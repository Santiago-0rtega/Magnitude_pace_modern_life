required <- c("metafor", "ggplot2", "dplyr", "here")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
parse(file = "remote_fit_one_small_study.R")
parse(file = "remote_finalize_small_study.R")
cat("PREFLIGHT_OK\n")
