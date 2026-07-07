save_plot_dual <- function(plot,
                           filename,
                           width  = 180,
                           height = 120,
                           units  = "mm",
                           dpi    = 300) {
  pdf_dir <- here::here("Rdata", "figures", "pdf")
  png_dir <- here::here("Rdata", "figures", "png")

  if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
  if (!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE)

  ggplot2::ggsave(
    filename = file.path(pdf_dir, paste0(filename, ".pdf")),
    plot   = plot,
    width  = width,
    height = height,
    units  = units
  )

  ggplot2::ggsave(
    filename = file.path(png_dir, paste0(filename, ".png")),
    plot   = plot,
    width  = width,
    height = height,
    units  = units,
    dpi    = dpi
  )

  invisible(NULL)
}
