# Combine per-model sensitivity parts into the CSVs the book's ch07 reads.
setwd(here::here())
source(here::here("Scripts", "01_paths.R"))
library(here)
suppressMessages({ library(dplyr); library(readr); library(purrr) })

part_dir  <- dir_out("tables", "sensitivity", "parts")
table_dir <- dir_out("tables", "sensitivity")

read_parts <- function(suffix) {
  fs <- list.files(part_dir, pattern = paste0(suffix, "\\.csv$"), full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  purrr::map_dfr(fs, ~ readr::read_csv(.x, show_col_types = FALSE))
}

variants <- c(sens1_n40 = "sens1_n40", sens2_nophylo = "sens2_nophylo", sens3_sysid = "sens3_sysid")

for (v in variants) {
  fe   <- read_parts("_fe")
  diag <- read_parts("_diag")
  if (is.null(fe)) next
  fe_v <- fe |> dplyr::filter(grepl(paste0("^", v, "_"), model_id))
  if (nrow(fe_v) == 0) next
  readr::write_csv(dplyr::filter(fe_v, submodel == "location"),
                   file.path(table_dir, paste0(v, "_location.csv")))
  readr::write_csv(dplyr::filter(fe_v, submodel == "scale"),
                   file.path(table_dir, paste0(v, "_scale.csv")))
  if (!is.null(diag)) {
    diag_v <- diag |> dplyr::filter(grepl(paste0("^", v, "_"), model_id))
    if (nrow(diag_v) > 0)
      readr::write_csv(diag_v, file.path(table_dir, paste0(v, "_diagnostics.csv")))
  }
  message("aggregated ", v, ": ", length(unique(fe_v$model_id)), " models")
}
message("Combined CSVs written to ", table_dir)
print(list.files(table_dir, pattern = "\\.csv$"))
