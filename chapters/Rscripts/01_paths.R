dir_data        <- here::here("data")
dir_data_clean  <- here::here("Rdata", "data_clean")
dir_effect_sizes <- here::here("Rdata", "effect_sizes")
dir_phylogeny   <- here::here("Rdata", "phylogeny")
dir_models      <- here::here("Rdata", "models")
dir_models_sens <- here::here("Rdata", "models", "sensitivity")
dir_diagnostics <- here::here("Rdata", "diagnostics")
dir_tables      <- here::here("Rdata", "tables")
dir_tables_sens <- here::here("Rdata", "tables", "sensitivity")
dir_fig_pdf     <- here::here("Rdata", "figures", "pdf")
dir_fig_png     <- here::here("Rdata", "figures", "png")
dir_fig_jpg     <- here::here("Rdata", "figures", "jpg")

dirs_all <- c(
  dir_data_clean, dir_effect_sizes, dir_phylogeny,
  dir_models, dir_models_sens, dir_diagnostics,
  dir_tables, dir_tables_sens,
  dir_fig_pdf, dir_fig_png, dir_fig_jpg
)

for (d in dirs_all) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
