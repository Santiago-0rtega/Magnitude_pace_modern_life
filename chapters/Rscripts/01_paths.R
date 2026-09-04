# Single source of truth for the derived-output tree. Every script builds its
# output paths through dir_out(), so the directory name is changed in one place
# instead of in a hundred literal here::here("Rdata", ...) calls. Override with
# options(pace.out_root = "somewhere_else") before sourcing this file.
dir_out <- function(...) here::here(getOption("pace.out_root", "Rdata"), ...)

dir_data        <- here::here("data")
dir_data_clean  <- dir_out("data_clean")
dir_effect_sizes <- dir_out("effect_sizes")
dir_phylogeny   <- dir_out("phylogeny")
dir_models      <- dir_out("models")
dir_models_sens <- dir_out("models", "sensitivity")
dir_diagnostics <- dir_out("diagnostics")
dir_tables      <- dir_out("tables")
dir_tables_sens <- dir_out("tables", "sensitivity")
dir_fig_pdf     <- dir_out("figures", "pdf")
dir_fig_png     <- dir_out("figures", "png")
dir_fig_jpg     <- dir_out("figures", "jpg")

dirs_all <- c(
  dir_data_clean, dir_effect_sizes, dir_phylogeny,
  dir_models, dir_models_sens, dir_diagnostics,
  dir_tables, dir_tables_sens,
  dir_fig_pdf, dir_fig_png, dir_fig_jpg
)

for (d in dirs_all) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
