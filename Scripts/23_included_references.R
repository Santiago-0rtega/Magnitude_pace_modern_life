# List of source references included in the primary analyses (ESM section S8).
#
# The analysed dataset is rebuilt exactly as the model scripts do it
# (run_m00_intercept_only.R / remote_fit_one_primary.R): cached SAFE lnM
# estimates -> apply_phylo_name_map() -> prepare_phylo_and_data(). Its ref_id
# set is cross-checked against the fitted m00 brmsfit when that file exists.
# Citation text (`reference`) and `doi` come from the raw PROCEED v6.2 CSV.
#
# Outputs
#   Rdata/tables/included_references.csv
#   <Proceedings B submission folder>/included_references.tex  (LaTeX fragment)
#
# Run from the repository root:
#   Rscript Scripts/23_included_references.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
})

options(pace.read_only = TRUE)

source(here::here("Scripts", "01_paths.R"))
source(here::here("Scripts", "02_read_clean_data.R"))
source(here::here("Scripts", "05_phylogeny.R"))

# Optional LaTeX fragment for a local manuscript folder: set PACE_MANUSCRIPT_DIR to enable.
tex_path <- file.path(
  Sys.getenv("PACE_MANUSCRIPT_DIR", unset = here::here("manuscript")),
  "included_references.tex"
)

# ── 1. Analysed dataset (primary model m00) ──────────────────────────────────
dat_es   <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
name_map <- readRDS(dir_out("phylogeny", "proceed_name_map.rds"))
A_full   <- readRDS(dir_out("phylogeny", "proceed_A_matrix.rds"))
dat_es   <- apply_phylo_name_map(dat_es, name_map)
dat_mod  <- suppressMessages(
  prepare_phylo_and_data(dat_es, A_full, fill_unmatched = TRUE, label = "m00")
)$dat_model

ref_ids <- sort(unique(as.character(dat_mod$ref_id)))
k_contrasts <- nrow(dat_mod)

m00_path <- dir_out("models", "m00_ls_intercept_only.rds")
if (file.exists(m00_path)) {
  m00_data <- readRDS(m00_path)$data
  if (nrow(m00_data) != k_contrasts ||
      !setequal(as.character(m00_data$ref_id), ref_ids)) {
    warning("Analysed dataset differs from the data stored in the m00 fit (",
            nrow(m00_data), " rows, ",
            dplyr::n_distinct(m00_data$ref_id), " references).")
  } else {
    message("Cross-check passed: m00 fit data has ", nrow(m00_data),
            " contrasts and the same ", length(ref_ids), " references.")
  }
}

if (length(ref_ids) != 254 || k_contrasts != 7186) {
  warning("Counts differ from manuscript: ", length(ref_ids),
          " references (stated 254), ", k_contrasts,
          " contrasts (stated 7,186).")
}

# ── 2. Citation text and DOI from PROCEED ─────────────────────────────────────
refs_raw <- read_proceed() |>
  dplyr::filter(ref_id %in% ref_ids) |>
  dplyr::distinct(ref_id, reference, doi)

if (anyDuplicated(refs_raw$ref_id)) {
  warning("Some ref_id values carry more than one reference/DOI string; ",
          "the first occurrence is used.")
  refs_raw <- dplyr::distinct(refs_raw, ref_id, .keep_all = TRUE)
}
stopifnot(setequal(refs_raw$ref_id, ref_ids))

# ── 3. Text repair and LaTeX escaping ─────────────────────────────────────────
# The PROCEED CSV contains literal U+FFFD replacement characters where the
# original characters were lost upstream. Repairs below were checked against
# Crossref author records where a DOI exists (Kristjansson, Moller, Lo Cascio
# Saetre); digit-FFFD-digit is a page-range dash; the Oke 2020 entry used FFFD
# as a field separator.
fffd <- "\uFFFD"
repair_text <- function(x) {
  x <- str_replace_all(x, paste0("Kristj", fffd, "nsson"), "Kristj\u00e1nsson")
  x <- str_replace_all(x, paste0("R", fffd, "s", fffd, "nen"), "R\u00e4s\u00e4nen")
  x <- str_replace_all(x, paste0("M", fffd, "ller 2006"), "M\u00f8ller 2006")
  x <- str_replace_all(x, paste0("S", fffd, "tre"), "S\u00e6tre")
  x <- str_replace_all(x, paste0("(\\d)", fffd, "(\\d)"), "\\1\u2013\\2")
  x <- str_replace_all(x, paste0("(\\d)", fffd, "(\\D)"), "\\1, \\2")
  x <- str_replace_all(x, paste0("(\\D)", fffd, "(\\d)"), "\\1, \\2")
  x
}

clean_ws <- function(x) str_squish(str_replace_all(x, "[\r\n\t\u00a0]+", " "))

latex_escape <- function(x) {
  x <- str_replace_all(x, fixed("\\"), "\u0001")
  x <- str_replace_all(x, "([&%$#_{}])", "\\\\\\1")
  x <- str_replace_all(x, fixed("~"), "\\textasciitilde{}")
  x <- str_replace_all(x, fixed("^"), "\\textasciicircum{}")
  x <- str_replace_all(x, fixed("\u0001"), "\\textbackslash{}")
  x
}

# Characters pdfLaTeX (utf8 inputenc, OT1/T1) may not handle, or that are
# better as ASCII ligatures. Common accented Latin letters stay as UTF-8.
char_map <- c(
  "\u2013" = "--", "\u2014" = "---", "\u2012" = "--", "\u2212" = "--",
  "\u2018" = "`",  "\u2019" = "'",   "\u201c" = "``", "\u201d" = "''",
  "\u00f8" = "{\\o}", "\u00d8" = "{\\O}", "\u00e6" = "{\\ae}", "\u00c6" = "{\\AE}",
  "\u00df" = "{\\ss}", "\u2026" = "\\ldots{}"
)
safe_utf8 <- "[\u00c0-\u00d6\u00d9-\u00dd\u00e0-\u00f6\u00f9-\u00fd\u00ff]"

to_latex_chars <- function(x) {
  for (k in names(char_map)) x <- str_replace_all(x, fixed(k), char_map[[k]])
  x
}

refs <- refs_raw |>
  dplyr::mutate(
    reference_clean = clean_ws(repair_text(reference)),
    doi_clean = doi |>
      clean_ws() |>
      str_remove_all("(?i)https?://(dx\\.)?doi\\.org/") |>
      str_remove("(?i)^doi:\\s*") |>
      dplyr::na_if("")
  )

# Report characters changed or still non-ASCII.
nonascii_before <- unique(unlist(str_extract_all(
  paste(refs$reference, refs$doi), "[^\\x01-\\x7f]")))
replaced <- intersect(names(char_map), unique(unlist(str_extract_all(
  refs$reference_clean, "[^\\x01-\\x7f]"))))
left_utf8 <- setdiff(unique(unlist(str_extract_all(
  refs$reference_clean, "[^\\x01-\\x7f]"))), names(char_map))
bad_left <- left_utf8[!str_detect(left_utf8, safe_utf8)]
if (any(str_detect(refs$reference_clean, fffd))) {
  warning("Unrepaired U+FFFD characters remain in: ",
          paste(refs$ref_id[str_detect(refs$reference_clean, fffd)], collapse = ", "))
}
if (length(bad_left) > 0) {
  warning("Non-ASCII characters not covered by the LaTeX map: ",
          paste(bad_left, collapse = " "))
}

# Alphabetical by first author (accent- and case-insensitive).
sort_key <- refs$reference_clean |>
  iconv(from = "UTF-8", to = "ASCII//TRANSLIT", sub = "") |>
  str_to_lower() |>
  str_remove_all("[^a-z0-9 ]")
refs <- refs[order(sort_key, refs$ref_id), ]

fmt_doi <- function(d) {
  if (is.na(d)) return("")
  # Split only where a new DOI starts; SICI DOIs contain ";" internally.
  parts <- str_trim(str_split(d, ";\\s*(?=10\\.)")[[1]])
  parts <- parts[nzchar(parts)]
  paste0(" doi: ", paste0("\\url{https://doi.org/", parts, "}", collapse = "; "))
}

refs <- refs |>
  dplyr::mutate(
    text_tex = to_latex_chars(latex_escape(reference_clean)),
    text_tex = ifelse(str_detect(text_tex, "\\.$"), text_tex, paste0(text_tex, ".")),
    entry    = paste0("\\item ", text_tex, vapply(doi_clean, fmt_doi, character(1)))
  )

# ── 4. Write outputs ──────────────────────────────────────────────────────────
readr::write_csv(
  refs |> dplyr::transmute(ref_ID = ref_id, reference, DOI = doi,
                           reference_clean),
  dir_out("tables", "included_references.csv")
)

fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)
tex <- c(
  "% Generated by Scripts/23_included_references.R -- do not edit by hand.",
  "\\subsection{S8. Studies included in the analyses}",
  "",
  sprintf(paste0("The analyses included %s contrasts from the following %s ",
                 "source references compiled in PROCEED v6.2. ",
                 "References are reproduced as recorded in PROCEED, ",
                 "including any incomplete or inconsistent bibliographic details."),
          fmt_n(k_contrasts), fmt_n(nrow(refs))),
  "",
  "\\begin{list}{}{\\setlength{\\leftmargin}{1.5em}\\setlength{\\itemindent}{-1.5em}}",
  refs$entry,
  "\\end{list}",
  ""
)
if (dir.exists(dirname(tex_path))) {
  con <- file(tex_path, open = "w", encoding = "UTF-8")
  writeLines(tex, con, useBytes = FALSE)
  close(con)
  message("Wrote: ", tex_path)
} else {
  warning("Submission folder not found; LaTeX fragment not written: ", tex_path)
}

# ── 5. Report ─────────────────────────────────────────────────────────────────
coauthors <- c("Ortega", "Santos", "Gotanda", "Sanderson", "Gorn\u00e9", "Gorne",
               "Hendry", "Nakagawa")
co_pat <- paste0("\\b(", paste(coauthors, collapse = "|"), ")\\b")
dup_doi <- refs |>
  dplyr::filter(!is.na(doi_clean)) |>
  dplyr::group_by(doi_clean) |>
  dplyr::filter(dplyr::n() > 1) |>
  dplyr::ungroup()
odd <- refs |>
  dplyr::filter(!str_detect(reference_clean, "\\d+\\s*(\\(\\d+\\))?\\s*[:,]\\s*e?\\d+") |
                  nchar(reference_clean) < 25)

cat("\n=== Included references ===\n")
cat("Contrasts:", k_contrasts, "| references:", nrow(refs),
    "| species:", dplyr::n_distinct(as.character(dat_mod$sp_ncbi)), "\n")
cat("Missing DOI:", sum(is.na(refs$doi_clean)), "\n")
cat("Empty reference text:", sum(is.na(refs$reference_clean) |
                                   !nzchar(refs$reference_clean)), "\n")
cat("Non-ASCII characters in source:", paste(nonascii_before, collapse = " "), "\n")
cat("Replaced by LaTeX macros/ASCII:", paste(replaced, collapse = " "), "\n")
cat("Left as UTF-8:", paste(setdiff(left_utf8, replaced), collapse = " "), "\n")
cat("\nEntries with U+FFFD repaired:\n")
print(as.data.frame(refs |> dplyr::filter(str_detect(reference, fffd)) |>
                      dplyr::select(ref_id, reference_clean)))
cat("\nEntries containing co-author surnames:\n")
print(as.data.frame(refs |> dplyr::filter(str_detect(reference_clean, co_pat)) |>
                      dplyr::select(ref_id, reference_clean, doi_clean)))
cat("\nDOIs shared by more than one ref_id:\n")
print(as.data.frame(dup_doi |> dplyr::select(ref_id, reference_clean, doi_clean)))
cat("\nEntries without volume/page pattern (possibly incomplete):\n")
print(as.data.frame(odd |> dplyr::select(ref_id, reference_clean, doi_clean)))
cat("\nSample entries:\n")
writeLines(refs$entry[c(1, round(nrow(refs) / 2), nrow(refs))])
