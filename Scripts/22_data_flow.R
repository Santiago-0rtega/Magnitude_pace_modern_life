# Data-flow (PRISMA-style) record flow for the ESM.
#
# Counts contrasts, source references and species at each sequential step of
# the analysis pipeline, using the project's own functions so that counts match
# the cached analysis objects exactly:
#   read_proceed() / clean_proceed()        (02_read_clean_data.R)
#   require_lnm_vars()                      (03_filters.R)
#   filter_generations(max_gen = 300)       (03_filters.R)
#   compute_or_read_lnm_safe()  [cache]     (04_lnm_safe.R)
#   apply_phylo_name_map() +
#   prepare_phylo_and_data()                (05_phylogeny.R)
#
# Outputs
#   Rdata/tables/data_flow_counts.csv
#   Figures/data_flow.pdf
#   <Proceedings B submission folder>/figures/figS7.pdf (copy, if folder exists)
#
# Run from the repository root:
#   Rscript Scripts/22_data_flow.R
#
# Species counts: steps up to "eligible" count distinct PROCEED `sp_ncbi`
# names; the analysed step counts distinct canonical tree-tip names
# (`sp_ncbi` after apply_phylo_name_map()/prepare_phylo_and_data()), which is
# what the primary models use. Synonym resolution can merge PROCEED names, so
# the two bases are not directly subtractable.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# Never recompute SAFE lnM or rebuild the phylogeny here; read caches only.
options(pace.read_only = TRUE)

source(here::here("Scripts", "01_paths.R"))
source(here::here("Scripts", "02_read_clean_data.R"))
source(here::here("Scripts", "03_filters.R"))
source(here::here("Scripts", "04_lnm_safe.R"))
source(here::here("Scripts", "05_phylogeny.R"))

# Optional copy into a local manuscript folder: set PACE_MANUSCRIPT_DIR to enable.
esm_fig_path <- file.path(
  Sys.getenv("PACE_MANUSCRIPT_DIR", unset = here::here("manuscript")),
  "figures", "figS7.pdf"
)

# ── 1. Pipeline steps ─────────────────────────────────────────────────────────
dat_raw   <- read_proceed()
dat_clean <- clean_proceed(dat_raw)
dat_lnm   <- require_lnm_vars(dat_clean)
dat_gen   <- filter_generations(dat_lnm, max_gen = 300)
dat_input <- add_derived_vars(dat_gen)

# Cached SAFE lnM estimates (the function drops non-finite / non-positive
# estimates or variances). Read only; recompute = FALSE.
dat_es <- compute_or_read_lnm_safe(
  dat_input,
  out_path  = dir_out("effect_sizes", "proceed_lnm_safe.rds"),
  diag_path = dir_out("tables", "safe_lnm_diagnostics.csv"),
  recompute = FALSE
)

# Check the cache is a subset of the freshly filtered data.
stopifnot(all(dat_es$es_id %in% dat_input$es_id))
if (file.exists(dir_out("data_clean", "proceed_clean_filtered.rds"))) {
  cached_clean <- readRDS(dir_out("data_clean", "proceed_clean_filtered.rds"))
  if (!setequal(cached_clean$es_id, dat_input$es_id)) {
    warning("Freshly filtered data differ from the cached ",
            "proceed_clean_filtered.rds; counts below use the fresh filter.")
  }
}

name_map <- readRDS(dir_out("phylogeny", "proceed_name_map.rds"))
A_full   <- readRDS(dir_out("phylogeny", "proceed_A_matrix.rds"))
dat_map  <- apply_phylo_name_map(dat_es, name_map)
dat_mod  <- suppressMessages(
  prepare_phylo_and_data(dat_map, A_full, fill_unmatched = TRUE,
                         label = "m00")
)$dat_model

# ── 2. Counts ─────────────────────────────────────────────────────────────────
n_ref <- function(d) dplyr::n_distinct(as.character(d$ref_id))
n_sp  <- function(d) {
  s <- as.character(d$sp_ncbi)
  dplyr::n_distinct(s[!is.na(s) & nzchar(s)])
}
lost_refs <- function(before, after) {
  length(setdiff(unique(as.character(before$ref_id)),
                 unique(as.character(after$ref_id))))
}

# Breakdown of the generations filter (applied as one function).
gen_excl <- dat_lnm |> dplyr::filter(!es_id %in% dat_gen$es_id)
n_gen_na   <- sum(is.na(gen_excl$generations))
n_gen_le0  <- sum(!is.na(gen_excl$generations) & gen_excl$generations <= 0)
n_gen_gt   <- sum(!is.na(gen_excl$generations) & gen_excl$generations > 300)
stopifnot(n_gen_na + n_gen_le0 + n_gen_gt == nrow(gen_excl))

# SAFE exclusions and their summary-statistic signature.
safe_excl <- dat_input |> dplyr::filter(!es_id %in% dat_es$es_id)
n_safe_zero_sd <- sum(safe_excl$sd1 == 0 | safe_excl$sd2 == 0)
safe_reason <- if (nrow(safe_excl) > 0 && n_safe_zero_sd == nrow(safe_excl)) {
  "SAFE lnM not estimable (within-sample SD of 0)"
} else {
  "SAFE lnM estimate or variance not finite or not positive"
}

# Phylogeny exclusions.
phy_excl <- dat_map |>
  dplyr::filter(is.na(sp_ncbi_canonical) | !nzchar(sp_ncbi_canonical))
n_phy_taxa <- dplyr::n_distinct(as.character(phy_excl$sp_ncbi))
n_phy_refs <- n_ref(phy_excl)

counts <- tibble::tribble(
  ~step_order, ~box, ~step, ~contrasts, ~references, ~species, ~species_basis, ~references_removed_entirely, ~note,
  1L, "main", "PROCEED v6.2 rates database (all records)",
    nrow(dat_raw), n_ref(dat_raw), n_sp(dat_raw), "PROCEED sp_ncbi names", NA_integer_, "",
  2L, "excluded", "Mean, SD or n missing for one or both samples (require_lnm_vars)",
    nrow(dat_clean) - nrow(dat_lnm), n_ref(dplyr::filter(dat_clean, !es_id %in% dat_lnm$es_id)), NA_integer_, "", lost_refs(dat_clean, dat_lnm), "",
  3L, "main", "Mean, SD and n reported for both samples",
    nrow(dat_lnm), n_ref(dat_lnm), n_sp(dat_lnm), "PROCEED sp_ncbi names", NA_integer_, "",
  4L, "excluded", "Elapsed generations unknown, not positive, or > 300 (filter_generations)",
    nrow(gen_excl), n_ref(gen_excl), NA_integer_, "", lost_refs(dat_lnm, dat_gen),
    sprintf("unknown = %d; <= 0 = %d; > 300 = %d", n_gen_na, n_gen_le0, n_gen_gt),
  5L, "main", "Known, positive elapsed generations <= 300",
    nrow(dat_gen), n_ref(dat_gen), n_sp(dat_gen), "PROCEED sp_ncbi names", NA_integer_, "",
  6L, "excluded", safe_reason,
    nrow(safe_excl), n_ref(safe_excl), NA_integer_, "", lost_refs(dat_input, dat_es),
    sprintf("zero SD in at least one sample = %d; n = 1 in a sample = %d",
            n_safe_zero_sd, sum(safe_excl$n1 == 1 | safe_excl$n2 == 1)),
  7L, "main", "Eligible contrasts (SAFE lnM estimated)",
    nrow(dat_es), n_ref(dat_es), n_sp(dat_es), "PROCEED sp_ncbi names", NA_integer_, "",
  8L, "excluded", "Species name not resolved to the phylogeny",
    nrow(phy_excl), n_phy_refs, n_phy_taxa, "PROCEED sp_ncbi names", lost_refs(dat_es, dat_mod),
    paste(sort(unique(as.character(phy_excl$sp_ncbi))), collapse = "; "),
  9L, "main", "Analysed in primary models (m00)",
    nrow(dat_mod), n_ref(dat_mod), n_sp(dat_mod), "canonical tree-tip names", NA_integer_,
    sprintf("min n1 + n2 = %s", format(min(dat_mod$n1 + dat_mod$n2)))
)

readr::write_csv(counts, dir_out("tables", "data_flow_counts.csv"))
print(as.data.frame(counts[, 1:8]))

# Reconcile against the numbers stated in the manuscript (report, never force).
expected <- c(eligible = 7237, excluded_phylo = 51, analysed = 7186,
              refs = 254, species = 253, min_n = 5)
observed <- c(eligible = nrow(dat_es), excluded_phylo = nrow(phy_excl),
              analysed = nrow(dat_mod), refs = n_ref(dat_mod),
              species = n_sp(dat_mod), min_n = min(dat_mod$n1 + dat_mod$n2))
mismatch <- observed != expected
if (any(mismatch)) {
  warning("Counts differ from manuscript: ",
          paste0(names(observed)[mismatch], " observed ", observed[mismatch],
                 " vs stated ", expected[mismatch], collapse = "; "))
} else {
  message("All manuscript counts reproduced: ",
          paste(names(observed), observed, sep = " = ", collapse = ", "))
}

# ── 3. Figure ─────────────────────────────────────────────────────────────────
fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
c_line <- function(d, sp_label = "species") {
  sprintf("%s contrasts  |  %s references  |  %s %s",
          fmt(nrow(d)), fmt(n_ref(d)), fmt(n_sp(d)), sp_label)
}

main_boxes <- tibble::tibble(
  id = 1:5,
  y  = c(9.0, 6.95, 4.9, 2.85, 0.8),
  label = c(
    paste0("PROCEED v6.2 rates database\n", c_line(dat_raw)),
    paste0("Mean, SD and n reported for both samples\n", c_line(dat_lnm)),
    paste0("Known, positive elapsed generations \u2264 300\n", c_line(dat_gen)),
    paste0("Eligible contrasts (SAFE lnM estimated)\n", c_line(dat_es)),
    paste0("Analysed in primary models\n", c_line(dat_mod))
  ),
  fontface = c("plain", "plain", "plain", "bold", "bold")
)

gen_lines <- c(
  if (n_gen_na  > 0) sprintf("generations unknown: %s", fmt(n_gen_na)),
  if (n_gen_le0 > 0) sprintf("generations \u2264 0: %s", fmt(n_gen_le0)),
  if (n_gen_gt  > 0) sprintf("generations > 300: %s", fmt(n_gen_gt))
)

excl_boxes <- tibble::tibble(
  id = 1:4,
  y  = (main_boxes$y[-5] + main_boxes$y[-1]) / 2,
  label = c(
    sprintf("Excluded: %s contrasts\nmean, SD or n missing for one\nor both samples",
            fmt(nrow(dat_clean) - nrow(dat_lnm))),
    paste0(sprintf("Excluded: %s contrasts\n", fmt(nrow(gen_excl))),
           paste(gen_lines, collapse = "\n")),
    sprintf("Excluded: %s contrasts\n%s",
            fmt(nrow(safe_excl)),
            if (grepl("SD of 0", safe_reason)) {
              "SAFE lnM not estimable\n(within-sample SD of 0)"
            } else {
              "SAFE lnM not finite or\nvariance not positive"
            }),
    sprintf("Excluded: %s contrasts\nspecies not resolved in the phylogeny\n(%s taxa, %s references)",
            fmt(nrow(phy_excl)), fmt(n_phy_taxa), fmt(n_phy_refs))
  )
)

main_x <- 3.3; main_w <- 5.8; main_h <- 0.95
excl_x <- 8.4; excl_w <- 3.8; excl_h <- 1.15

main_rect <- main_boxes |>
  dplyr::mutate(xmin = main_x - main_w / 2, xmax = main_x + main_w / 2,
                ymin = y - main_h / 2,      ymax = y + main_h / 2)
excl_rect <- excl_boxes |>
  dplyr::mutate(xmin = excl_x - excl_w / 2, xmax = excl_x + excl_w / 2,
                ymin = y - excl_h / 2,      ymax = y + excl_h / 2)

down_arrows <- tibble::tibble(
  x = main_x, xend = main_x,
  y = main_rect$ymin[-5], yend = main_rect$ymax[-1]
)
side_arrows <- tibble::tibble(
  x = main_x, xend = excl_rect$xmin,
  y = excl_boxes$y, yend = excl_boxes$y
)

arr <- grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")
txt_size <- 3.4

p_flow <- ggplot() +
  geom_segment(data = down_arrows,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arr, linewidth = 0.45, colour = "grey15") +
  geom_segment(data = side_arrows,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arr, linewidth = 0.45, colour = "grey15") +
  geom_rect(data = main_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", colour = "black", linewidth = 0.5) +
  geom_rect(data = excl_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey93", colour = "grey35", linewidth = 0.4) +
  geom_text(data = main_rect,
            aes(x = main_x, y = y, label = label, fontface = fontface),
            size = txt_size, lineheight = 1.15, colour = "black") +
  geom_text(data = excl_rect,
            aes(x = excl_x, y = y, label = label),
            size = txt_size * 0.95, lineheight = 1.1, colour = "grey10") +
  coord_cartesian(xlim = c(0.3, 10.4), ylim = c(0.15, 9.65), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

fig_path <- here::here("Figures", "data_flow.pdf")
dev_fun  <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
ggsave(fig_path, p_flow, width = 8, height = 7, units = "in", device = dev_fun)
message("Saved figure: ", fig_path)

if (dir.exists(dirname(esm_fig_path))) {
  ok <- file.copy(fig_path, esm_fig_path, overwrite = TRUE)
  message(if (ok) "Copied to: " else "Copy FAILED: ", esm_fig_path)
} else {
  warning("ESM figure folder not found; copy skipped: ", dirname(esm_fig_path))
}
