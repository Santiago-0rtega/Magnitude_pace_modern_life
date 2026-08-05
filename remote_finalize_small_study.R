setwd("/home/ortegara/Documents/PACE")
library(here)
library(ggplot2)
library(metafor)
library(patchwork)
source(here::here("R", "20_small_study_effects.R"))

out_dir <- here::here("outputs", "small_study")
se <- readRDS(file.path(out_dir, "n_se.rds"))
vv <- readRDS(file.path(out_dir, "n_v.rds"))
stopifnot(se$n == vv$n)

table_dir <- here::here("outputs", "tables")
model_dir <- here::here("outputs", "models")
figure_dir <- here::here("outputs", "figures", "publication")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

results <- rbind(
  read.csv(file.path(out_dir, "n_se.csv"), check.names = FALSE),
  read.csv(file.path(out_dir, "n_v.csv"), check.names = FALSE)
)
write.csv(results, file.path(table_dir, "small_study_effects.csv"), row.names = FALSE)
saveRDS(list(n = se$n, fit_se = se$fit, fit_v = vv$fit),
        file.path(model_dir, "small_study_effects.rds"))

inputs <- prepare_small_study_data(
  effect_size_path = here::here("outputs", "effect_sizes", "proceed_lnm_safe.rds"),
  A_path = here::here("outputs", "phylogeny", "proceed_A_matrix.rds"),
  name_map_path = here::here("outputs", "phylogeny", "proceed_name_map.rds")
)
p_se <- build_small_study_figure(
  se$fit, inputs$data, predictor = "n_se",
  bubble_fill = "#DC143C", line_colour = "#DC143C"
)
p_v <- build_small_study_figure(
  vv$fit, inputs$data, predictor = "n_v",
  bubble_fill = "#FF8C00", line_colour = "#8B4500"
)
p_se <- p_se + ggplot2::theme(legend.position = "bottom")
p_v <- p_v + ggplot2::theme(legend.position = "none")
p <- ((p_se + p_v) / patchwork::guide_area()) +
  patchwork::plot_layout(heights = c(1, 0.10), guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "center",
    plot.tag = ggplot2::element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.99)
  )

ggsave(file.path(figure_dir, "small_study_effects.png"), p,
       width = 13, height = 6, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "small_study_effects.pdf"), p,
       width = 13, height = 6, device = grDevices::cairo_pdf)
cat("FINALIZED:", se$n, "contrasts\n")
