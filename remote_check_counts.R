setwd(here::here())
source(here::here("chapters", "Rscripts", "01_paths.R"))
files <- c(
  list.files(dir_out("models"), pattern = "^m(00|01|02|03|04|05|07|11).*rds$", full.names = TRUE),
  list.files(dir_out("models", "primary"), pattern = "years|generations", full.names = TRUE)
)
for (file in files) {
  object <- readRDS(file)
  data <- object[["data"]]
  cat(file, nrow(data),
      sum(!is.na(data[["log10_years"]])),
      sum(!is.na(data[["log10_generations"]])), "\n")
}
effect_data <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
cat("effect-size data", nrow(effect_data),
    sum(!is.na(effect_data[["log10_years"]])),
    sum(!is.na(effect_data[["log10_generations"]])), "\n")
