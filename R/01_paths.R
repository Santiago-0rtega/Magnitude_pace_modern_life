dir_data        <- here::here("data")
dir_data_clean  <- here::here("outputs", "data_clean")
dir_effect_sizes <- here::here("outputs", "effect_sizes")
dir_phylogeny   <- here::here("outputs", "phylogeny")
dir_models      <- here::here("outputs", "models")
dir_models_sens <- here::here("outputs", "models", "sensitivity")
dir_diagnostics <- here::here("outputs", "diagnostics")
dir_tables      <- here::here("outputs", "tables")
dir_tables_sens <- here::here("outputs", "tables", "sensitivity")
dir_fig_pdf     <- here::here("outputs", "figures", "pdf")
dir_fig_png     <- here::here("outputs", "figures", "png")
dir_fig_jpg     <- here::here("outputs", "figures", "jpg")

dirs_all <- c(
  dir_data_clean, dir_effect_sizes, dir_phylogeny,
  dir_models, dir_models_sens, dir_diagnostics,
  dir_tables, dir_tables_sens,
  dir_fig_pdf, dir_fig_png, dir_fig_jpg
)

for (d in dirs_all) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
