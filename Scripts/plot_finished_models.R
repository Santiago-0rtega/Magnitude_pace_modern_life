setwd(here::here())
source(here::here("Scripts", "01_paths.R"))
library(here)
source(here::here("Scripts", "00_packages_fit.R"))
source(here::here("Scripts", "06_model_registry.R"))
source(here::here("Scripts", "09_model_summaries.R"))
source(here::here("Scripts", "11_plotting_orchard_like.R"))

out_dir <- dir_out("figures", "model_grid")
out_dir_pdf <- here::here("Figures", "model_grid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_pdf, recursive = TRUE, showWarnings = FALSE)

finished <- c("m01_ls_disturbance", "m02_ls_design", "m03_ls_log10_years",
              "m04_ls_log10_generations", "m05_ls_trait_type")

moderators  <- c("disturbance", "design", "log10_years", "log10_generations", "trait_type")
labels      <- c("Disturbance context", "Comparison design",
                 "Elapsed time (log10 years)", "Elapsed time (log10 generations)", "Trait type")
mod_types   <- c("categorical", "categorical", "continuous", "continuous", "categorical")

dat_es <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))

save_plot <- function(p, stem, width, height) {
  cairo_pdf(file.path(out_dir_pdf, paste0(stem, ".pdf")), width = width, height = height)
  print(p); dev.off()
  png(file.path(out_dir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 300, type = "cairo")
  print(p); dev.off()
  message("Saved: ", stem)
}

for (i in seq_along(finished)) {
  nm        <- finished[i]
  mid       <- sub("_ls_.*", "", nm)
  moderator <- moderators[i]
  label     <- labels[i]
  mod_type  <- mod_types[i]

  message("\nPlotting ", nm)
  fit      <- readRDS(dir_out("models", paste0(nm, ".rds")))
  fe       <- extract_fixed_effects(fit, mid, moderator)
  ls_split <- split_location_scale(fe)

  if (mod_type == "categorical") {
    p_loc <- plot_categorical_location(ls_split$location, label)
    p_scl <- plot_categorical_scale(ls_split$scale, label)
  } else {
    dat_model <- dat_es[!is.na(dat_es[[moderator]]), ]
    p_loc <- plot_continuous_location(fit, moderator, label, dat_model)
    p_scl <- plot_continuous_scale(fit, moderator, label, dat_model)
  }

  p_comb <- plot_combined_ls(p_loc, p_scl)
  h <- if (mod_type == "categorical") 7 else 6

  save_plot(p_loc,  paste0(mid, "_location"), width = 8, height = h * 0.55)
  save_plot(p_scl,  paste0(mid, "_scale"),    width = 8, height = h * 0.55)
  save_plot(p_comb, paste0(mid, "_combined"), width = 8, height = h)
}

message("\nAll plots saved to: ", out_dir)
