setwd(here::here())
library(here)
for (f in c("00_packages_fit", "01_paths", "02_read_clean_data", "03_filters", "04_lnm_safe"))
  source(here::here("Scripts", paste0(f, ".R")))

clean_path <- dir_out("data_clean", "proceed_clean_filtered.rds")
es_path <- dir_out("effect_sizes", "proceed_lnm_safe.rds")
diag_path <- dir_out("tables", "safe_lnm_diagnostics.csv")
for (p in c(clean_path, es_path, diag_path))
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

dat <- read_proceed() |> clean_proceed() |> require_lnm_vars() |>
  filter_generations(max_gen = 300) |> add_derived_vars()
stopifnot(nrow(dat) > 0, all(!is.na(dat$generations)),
          all(dat$generations > 0), all(dat$generations <= 300))
saveRDS(dat, clean_path)
cat("Filtered contrasts:", nrow(dat), "\n")

dat_es <- compute_or_read_lnm_safe(dat, out_path = es_path,
  diag_path = diag_path, recompute = TRUE, seed = 123)
stopifnot(all(!is.na(dat_es$generations)), all(dat_es$generations <= 300))
cat("SAFE effect sizes retained:", nrow(dat_es), "\n")
