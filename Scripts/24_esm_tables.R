# Electronic supplementary material tables S2-S4 as standalone LaTeX fragments.
#
#   Table S2  Marginal means of the primary location-scale models (m00, m01,
#             m02, m03, m04, m05, m07) with per-level k, references, species.
#   Table S3  Pairwise contrasts among moderator levels (m01, m02, m05, m07).
#   Table S4  Time-adjusted (additive, non-preregistered) models: elapsed-time
#             slopes and pairwise moderator contrasts, with convergence summary.
#
# Inputs (tracked): Rdata/tables/contrasts/m0X_{location,scale}_{emmeans,contrasts}.csv,
#   Rdata/tables/primary/primary_fixef.csv, Rdata/tables/additive_contrasts.csv
#   (from Scripts/remote_additive_contrasts.R), Rdata/tables/additive_diagnostics.csv.
# Inputs (local caches): effect sizes, phylogeny name map and A matrix; the
#   intercept-only fit (m00) for the overall estimates, cached afterwards in
#   Rdata/tables/esm_m00_fixef.csv. Fitted primary models, when present, are
#   used to cross-check the reconstructed analysed data (skip with
#   PACE_ESM_CHECK_FITS=0).
#
# Counts per model are rebuilt exactly as in Scripts/remote_fit_one_primary.R:
#   SAFE lnM cache -> apply_phylo_name_map() -> rows with a non-missing
#   moderator -> prepare_phylo_and_data().
#
# Contrast direction (S3 and S4): levels are listed in factor (reference-coding)
#   order; each pair is reported as later-listed level minus earlier-listed
#   level. S3 contrasts stored in the opposite direction are negated (the
#   equal-tailed interval bounds are negated and swapped, which is exact).
#
# Outputs
#   Rdata/tables/esm_table_S2.csv, esm_table_S3.csv, esm_table_S4.csv
#   <PACE_MANUSCRIPT_DIR>/esm_table_S2.tex, esm_table_S3.tex, esm_table_S4.tex
#
# Run from the repository root:
#   Rscript Scripts/24_esm_tables.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
})

options(pace.read_only = TRUE)
source(here::here("Scripts", "01_paths.R"))
source(here::here("Scripts", "05_phylogeny.R"))

# Optional LaTeX fragments for a local manuscript folder: set PACE_MANUSCRIPT_DIR.
tex_dir <- Sys.getenv("PACE_MANUSCRIPT_DIR", unset = here::here("manuscript"))
check_fits <- Sys.getenv("PACE_ESM_CHECK_FITS", unset = "1") != "0"

# ── 0. Labels ─────────────────────────────────────────────────────────────────
level_labels <- list(
  disturbance = c(
    "Climate change"            = "Climate change",
    "Hunt_harv"                 = "Hunting or harvesting",
    "Introduction"              = "Introduced populations",
    "Landscape change"          = "Landscape change",
    "Other"                     = "\\textit{In situ} natural variation",
    "Pollution"                 = "Pollution",
    "Response to introductions" = "Responses to introduced species"),
  design = c("Allochronic" = "Allochronic", "Synchronic" = "Synchronic"),
  trait_type = c(
    "behaviour"       = "Behaviour",
    "growth"          = "Growth",
    "otherLH"         = "Other life-history",
    "othermorphology" = "Other morphology",
    "phenology"       = "Phenology",
    "physio"          = "Physiology",
    "response"        = "Response",
    "size"            = "Body size"),
  genphen = c("Genetic" = "Common garden/QG", "Phenotypic" = "Field")
)
moderator_labels <- c(
  disturbance       = "Disturbance context",
  design            = "Comparison design",
  log10_years       = "Elapsed time ($\\log_{10}$ years)",
  log10_generations = "Elapsed time ($\\log_{10}$ generations)",
  trait_type        = "Trait type",
  genphen           = "Study type"
)
timevar_labels <- c(log10_years = "$\\log_{10}$ years",
                    log10_generations = "$\\log_{10}$ generations")

cat_models <- tibble::tribble(
  ~model_id, ~moderator,
  "m01", "disturbance",
  "m02", "design",
  "m05", "trait_type",
  "m07", "genphen"
)
cont_models <- tibble::tribble(
  ~model_id, ~moderator,         ~fixef_id,
  "m03", "log10_years",          "primary_log10_years",
  "m04", "log10_generations",    "primary_log10_generations"
)
model_files <- c(
  m00 = "m00_ls_intercept_only", m01 = "m01_ls_disturbance",
  m02 = "m02_ls_design", m03 = "m03_ls_log10_years",
  m04 = "m04_ls_log10_generations", m05 = "m05_ls_trait_type",
  m07 = "m07_ls_genphen"
)

# ── 1. Analysed data per model ────────────────────────────────────────────────
dat_es   <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
name_map <- readRDS(dir_out("phylogeny", "proceed_name_map.rds"))
A_full   <- readRDS(dir_out("phylogeny", "proceed_A_matrix.rds"))
dat_es   <- apply_phylo_name_map(dat_es, name_map)

model_data <- function(model_id, moderator = NULL, categorical = TRUE) {
  d <- dat_es
  if (!is.null(moderator)) {
    d <- d[!is.na(d[[moderator]]), ]
    if (categorical) d[[moderator]] <- droplevels(factor(d[[moderator]]))
  }
  d <- d |> dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
  suppressMessages(prepare_phylo_and_data(d, A_full, label = model_id))$dat_model
}

count_rows <- function(d) {
  c(k = nrow(d),
    n_refs = dplyr::n_distinct(as.character(d$ref_id)),
    n_species = dplyr::n_distinct(as.character(d$sp_ncbi)))
}

all_ids <- c("m00", cat_models$model_id, cont_models$model_id)
mods    <- c(m00 = NA, setNames(cat_models$moderator, cat_models$model_id),
             setNames(cont_models$moderator, cont_models$model_id))
dat_by_model <- lapply(all_ids, function(id) {
  mod <- mods[[id]]
  if (is.na(mod)) model_data(id)
  else model_data(id, mod, categorical = id %in% cat_models$model_id)
})
names(dat_by_model) <- all_ids

m00_counts <- count_rows(dat_by_model$m00)
if (!identical(unname(m00_counts), c(7186L, 254L, 253L)))
  warning("m00 counts differ from the manuscript (7,186 / 254 / 253): ",
          paste(m00_counts, collapse = " / "))

# Optional cross-check against the data stored in the fitted models.
fit_check <- list()
m00_fit <- NULL
for (id in all_ids) {
  path <- dir_out("models", paste0(model_files[[id]], ".rds"))
  if (!file.exists(path) || (!check_fits && id != "m00")) next
  fit <- readRDS(path)
  fd  <- fit$data
  d   <- dat_by_model[[id]]
  ok  <- nrow(fd) == nrow(d) &&
    isTRUE(all.equal(fd$yi_lnM_safe, d$yi_lnM_safe)) &&
    identical(as.character(fd$ref_id), as.character(d$ref_id)) &&
    identical(as.character(fd$sp_ncbi), as.character(d$sp_ncbi))
  fit_check[[id]] <- data.frame(model_id = id, fit_k = nrow(fd),
                                rebuilt_k = nrow(d), identical_rows = ok)
  if (id == "m00") m00_fit <- fit
  rm(fit); invisible(gc())
}
if (length(fit_check)) {
  fit_check <- do.call(rbind, fit_check)
  cat("\n=== Rebuilt analysed data vs fitted-model data ===\n")
  print(fit_check, row.names = FALSE)
  if (!all(fit_check$identical_rows))
    warning("Rebuilt analysed data differ from the fitted-model data for: ",
            paste(fit_check$model_id[!fit_check$identical_rows], collapse = ", "))
}

# ── 2. Estimates ──────────────────────────────────────────────────────────────
tab_dir <- dir_out("tables")
rd <- function(...) readr::read_csv(file.path(tab_dir, ...), show_col_types = FALSE)

# m00 overall estimates (posterior mean, 95% equal-tailed CrI); cached as CSV
m00_cache <- file.path(tab_dir, "esm_m00_fixef.csv")
if (!is.null(m00_fit)) {
  dr <- posterior::as_draws_df(m00_fit)
  sm <- function(x) c(estimate = mean(x),
                      lower = unname(stats::quantile(x, 0.025)),
                      upper = unname(stats::quantile(x, 0.975)))
  m00_est <- rbind(
    data.frame(component = "location", term = "Intercept", t(sm(dr$b_Intercept))),
    data.frame(component = "scale", term = "sigma_Intercept",
               t(sm(dr$b_sigma_Intercept))))
  readr::write_csv(m00_est, m00_cache)
  rm(m00_fit, dr); invisible(gc())
} else if (file.exists(m00_cache)) {
  m00_est <- readr::read_csv(m00_cache, show_col_types = FALSE)
} else stop("Neither the m00 fit nor ", m00_cache, " is available.")

fixef <- rd("primary", "primary_fixef.csv")

# ── 3. Formatting helpers ─────────────────────────────────────────────────────
fmt_num <- function(x) {
  s <- sprintf("%.3f", x)
  s
}
fmt_ci <- function(est, lo, hi, mark = FALSE) {
  paste0("$", fmt_num(est), "\\ [", fmt_num(lo), ",\\ ", fmt_num(hi), "]",
         ifelse(mark, "^{\\dagger}", ""), "$")
}
fmt_int <- function(x) ifelse(x >= 1000,
                              sub("(\\d)(\\d{3})$", "\\1{,}\\2", as.character(x)),
                              as.character(x))
excl0 <- function(lo, hi) lo > 0 | hi < 0

lab <- function(moderator, code) {
  m <- level_labels[[moderator]]
  stopifnot(all(code %in% names(m)))
  unname(m[code])
}
plain <- function(x) {      # LaTeX label -> plain text for CSV/console
  x <- gsub("\\\\textit\\{([^}]*)\\}", "\\1", x)
  x <- gsub("\\$\\\\log_\\{10\\}\\$", "log10", x)
  x <- gsub("\\$-\\$", "-", x)
  x
}

ltx_begin <- function(colspec, caption, label, header, ncol) {
  c("\\begingroup",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\setlength{\\LTcapwidth}{\\textwidth}",
    paste0("\\begin{longtable}{", colspec, "}"),
    paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
    "\\toprule", header, "\\midrule",
    "\\endfirsthead",
    paste0("\\multicolumn{", ncol, "}{@{}l}{\\textit{\\tablename~\\thetable{} (continued)}}\\\\"),
    "\\toprule", header, "\\midrule",
    "\\endhead",
    "\\midrule",
    paste0("\\multicolumn{", ncol, "}{r@{}}{\\textit{Continued on next page}}\\\\"),
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot")
}
ltx_end <- c("\\end{longtable}", "\\endgroup")

write_tex <- function(lines, file, script = "Scripts/24_esm_tables.R") {
  path <- file.path(tex_dir, file)
  lines <- c(paste0("% Generated by ", script, " -- do not edit by hand."), lines, "")
  if (dir.exists(tex_dir)) {
    con <- file(path, open = "w", encoding = "UTF-8")
    writeLines(lines, con, useBytes = FALSE)
    close(con)
    message("Wrote: ", file)
  } else {
    warning("Manuscript folder not found; LaTeX fragment not written: ", file)
  }
  invisible(lines)
}

check_balance <- function(lines, name) {
  txt <- paste(lines, collapse = "\n")
  txt <- gsub("\\\\[{}]", "", txt)
  opens  <- lengths(regmatches(txt, gregexpr("\\{", txt)))
  closes <- lengths(regmatches(txt, gregexpr("\\}", txt)))
  envs_b <- regmatches(txt, gregexpr("\\\\begin\\{[^}]+\\}", txt))[[1]]
  envs_e <- regmatches(txt, gregexpr("\\\\end\\{[^}]+\\}", txt))[[1]]
  dollars <- lengths(regmatches(txt, gregexpr("\\$", txt)))
  grp <- c(lengths(regmatches(txt, gregexpr("\\\\begingroup", txt))),
           lengths(regmatches(txt, gregexpr("\\\\endgroup", txt))))
  ok <- opens == closes &&
    identical(sort(sub("begin", "", envs_b)), sort(sub("end", "", envs_e))) &&
    dollars %% 2 == 0 && grp[1] == grp[2]
  cat(sprintf("Balance %s: braces %d/%d, environments %d/%d, $ %d, groups %d/%d -> %s\n",
              name, opens, closes, length(envs_b), length(envs_e), dollars,
              grp[1], grp[2], if (ok) "OK" else "FAILED"))
  if (!ok) stop("LaTeX balance check failed for ", name)
}

orient_contrasts <- function(df, levels_order) {
  # Re-express "A - B" rows as later-listed minus earlier-listed level.
  parts <- strsplit(df$contrast, " - ", fixed = TRUE)
  a <- vapply(parts, `[`, "", 1); b <- vapply(parts, `[`, "", 2)
  stopifnot(all(c(a, b) %in% levels_order))
  flip <- match(a, levels_order) < match(b, levels_order)
  est <- ifelse(flip, -df$estimate, df$estimate)
  lo  <- ifelse(flip, -df$upper.CrI, df$lower.CrI)
  hi  <- ifelse(flip, -df$lower.CrI, df$upper.CrI)
  minuend <- ifelse(flip, b, a); subtrahend <- ifelse(flip, a, b)
  out <- data.frame(minuend, subtrahend, estimate = est, lower = lo, upper = hi,
                    flipped = flip, stringsAsFactors = FALSE)
  out[order(match(out$subtrahend, levels_order), match(out$minuend, levels_order)), ]
}

# ── 4. Table S2: marginal means ───────────────────────────────────────────────
s2_rows <- list()
add_s2 <- function(block, row_label, is_header, counts, loc = NULL, scl = NULL) {
  s2_rows[[length(s2_rows) + 1]] <<- data.frame(
    block = block, row = row_label, header = is_header,
    k = counts[["k"]], n_refs = counts[["n_refs"]], n_species = counts[["n_species"]],
    loc_est = if (is.null(loc)) NA else loc[[1]], loc_lo = if (is.null(loc)) NA else loc[[2]],
    loc_hi = if (is.null(loc)) NA else loc[[3]],
    scl_est = if (is.null(scl)) NA else scl[[1]], scl_lo = if (is.null(scl)) NA else scl[[2]],
    scl_hi = if (is.null(scl)) NA else scl[[3]], stringsAsFactors = FALSE)
}

# Overall (m00)
l0 <- m00_est[m00_est$component == "location", ]; s0 <- m00_est[m00_est$component == "scale", ]
add_s2("Overall (intercept-only model)", "Overall (intercept-only model)", TRUE, m00_counts)
add_s2("Overall (intercept-only model)", "Intercept", FALSE,
       c(k = NA, n_refs = NA, n_species = NA),
       c(l0$estimate, l0$lower, l0$upper), c(s0$estimate, s0$lower, s0$upper))

add_cat_block <- function(model_id, moderator) {
  d <- dat_by_model[[model_id]]
  levs <- levels(d[[moderator]])
  le <- rd("contrasts", paste0(model_id, "_location_emmeans.csv"))
  se <- rd("contrasts", paste0(model_id, "_scale_emmeans.csv"))
  stopifnot(setequal(le$level, levs), setequal(se$level, levs))
  blk <- moderator_labels[[moderator]]
  add_s2(blk, blk, TRUE, count_rows(d))
  for (lv in levs) {
    l <- le[le$level == lv, ]; s <- se[se$level == lv, ]
    add_s2(blk, lab(moderator, lv), FALSE, count_rows(d[d[[moderator]] == lv, ]),
           c(l$emmean, l$lower.CrI, l$upper.CrI), c(s$emmean, s$lower.CrI, s$upper.CrI))
  }
}
add_cont_block <- function(model_id, moderator, fixef_id) {
  d  <- dat_by_model[[model_id]]
  fe <- fixef[fixef$model_id == fixef_id, ]
  g  <- function(term) { r <- fe[fe$term == term, ]; stopifnot(nrow(r) == 1)
                         c(r$estimate, r$q2_5, r$q97_5) }
  blk <- moderator_labels[[moderator]]
  add_s2(blk, blk, TRUE, count_rows(d))
  na3 <- c(k = NA, n_refs = NA, n_species = NA)
  add_s2(blk, "Intercept", FALSE, na3, g("Intercept"), g("sigma_Intercept"))
  add_s2(blk, "Slope", FALSE, na3, g(moderator), g(paste0("sigma_", moderator)))
}

add_cat_block("m01", "disturbance")
add_cat_block("m02", "design")
add_cont_block("m03", "log10_years", "primary_log10_years")
add_cont_block("m04", "log10_generations", "primary_log10_generations")
add_cat_block("m05", "trait_type")
add_cat_block("m07", "genphen")
s2 <- do.call(rbind, s2_rows)

# Per-level counts must add up to the model totals (k only; references and
# species can be shared among levels).
for (b in unique(s2$block[s2$header & s2$block %in% moderator_labels[cat_models$moderator]])) {
  sb <- s2[s2$block == b, ]
  stopifnot(sum(sb$k[!sb$header]) == sb$k[sb$header])
}

fmt_cnt <- function(x) ifelse(is.na(x), "", fmt_int(x))
s2_tex_rows <- character(0)
for (i in seq_len(nrow(s2))) {
  r <- s2[i, ]
  if (r$header) {
    if (length(s2_tex_rows)) s2_tex_rows <- c(s2_tex_rows, "\\addlinespace")
    s2_tex_rows <- c(s2_tex_rows, paste0(
      "\\multicolumn{1}{@{}l}{\\textbf{", r$row, "}} & ", fmt_cnt(r$k), " & ",
      fmt_cnt(r$n_refs), " & ", fmt_cnt(r$n_species), " & & \\\\"))
  } else {
    s2_tex_rows <- c(s2_tex_rows, paste0(
      "\\quad ", r$row, " & ", fmt_cnt(r$k), " & ", fmt_cnt(r$n_refs), " & ",
      fmt_cnt(r$n_species), " & ", fmt_ci(r$loc_est, r$loc_lo, r$loc_hi), " & ",
      fmt_ci(r$scl_est, r$scl_lo, r$scl_hi), " \\\\"))
  }
}
s2_caption <- paste0(
  "\\textbf{Marginal means of the primary location--scale models.} ",
  "Values are posterior means with 95\\% equal-tailed credible intervals (CrI). ",
  "Location estimates are on the $\\ln M$ scale; scale estimates are the natural ",
  "logarithm of the residual standard deviation, $\\ln\\sigma$. ",
  "For categorical moderators, each row gives the population-level prediction ",
  "for that level (random effects set to zero). For the intercept-only and ",
  "elapsed-time models, rows give the intercept and the slope per one-unit ",
  "increase in the $\\log_{10}$-transformed elapsed time. ",
  "$k$, number of contrasts; Refs, number of source references; Spp, number ",
  "of species (accepted species-level names after taxonomic reconciliation). ",
  "Bold rows give these counts for all data analysed by the model; indented rows ",
  "give them for each level. References and species can contribute to more than ",
  "one level.")
s2_header <- paste0("Model / level & $k$ & Refs & Spp & ",
                    "Location, $\\ln M$ [95\\% CrI] & Scale, $\\ln\\sigma$ [95\\% CrI] \\\\")
s2_tex <- c(ltx_begin("@{}lrrrcc@{}", s2_caption, "tab:S2-marginal-means", s2_header, 6),
            s2_tex_rows, ltx_end)

# ── 5. Table S3: pairwise contrasts (primary models) ──────────────────────────
s3_rows <- list()
for (i in seq_len(nrow(cat_models))) {
  id <- cat_models$model_id[i]; mod <- cat_models$moderator[i]
  levs <- levels(dat_by_model[[id]][[mod]])
  lc <- orient_contrasts(rd("contrasts", paste0(id, "_location_contrasts.csv")), levs)
  sc <- orient_contrasts(rd("contrasts", paste0(id, "_scale_contrasts.csv")), levs)
  stopifnot(nrow(lc) == choose(length(levs), 2), nrow(sc) == nrow(lc))
  m <- merge(lc, sc, by = c("minuend", "subtrahend"), suffixes = c("_loc", "_scl"),
             sort = FALSE)
  m <- m[order(match(m$subtrahend, levs), match(m$minuend, levs)), ]
  s3_rows[[i]] <- data.frame(block = moderator_labels[[mod]], moderator = mod,
                             m, stringsAsFactors = FALSE)
}
s3 <- do.call(rbind, s3_rows)
s3$contrast_label <- paste0(mapply(lab, s3$moderator, s3$minuend),
                            " $-$ ", mapply(lab, s3$moderator, s3$subtrahend))

s3_tex_rows <- character(0)
for (b in unique(s3$block)) {
  if (length(s3_tex_rows)) s3_tex_rows <- c(s3_tex_rows, "\\addlinespace")
  s3_tex_rows <- c(s3_tex_rows, paste0("\\multicolumn{3}{@{}l}{\\textbf{", b, "}} \\\\"))
  sb <- s3[s3$block == b, ]
  s3_tex_rows <- c(s3_tex_rows, paste0(
    "\\quad ", sb$contrast_label, " & ",
    fmt_ci(sb$estimate_loc, sb$lower_loc, sb$upper_loc, excl0(sb$lower_loc, sb$upper_loc)),
    " & ",
    fmt_ci(sb$estimate_scl, sb$lower_scl, sb$upper_scl, excl0(sb$lower_scl, sb$upper_scl)),
    " \\\\"))
}
s3_caption <- paste0(
  "\\textbf{Pairwise contrasts among moderator levels in the primary location--scale models.} ",
  "Values are posterior means of the differences between level-specific estimates, ",
  "with 95\\% equal-tailed credible intervals (CrI). Each contrast is the first-named ",
  "level minus the second-named level. Location contrasts are differences in $\\ln M$; ",
  "scale contrasts are differences in the natural logarithm of the residual standard ",
  "deviation, $\\ln\\sigma$. $^{\\dagger}$95\\% CrI excludes zero.")
s3_header <- "Contrast & Location, $\\Delta\\ln M$ [95\\% CrI] & Scale, $\\Delta\\ln\\sigma$ [95\\% CrI] \\\\"
s3_tex <- c(ltx_begin("@{}lcc@{}", s3_caption, "tab:S3-contrasts", s3_header, 3),
            s3_tex_rows, ltx_end)

# ── 6. Table S4: time-adjusted (additive) models ──────────────────────────────
add  <- rd("additive_contrasts.csv")
diag <- rd("additive_diagnostics.csv")
nonconv <- "trait_type_plus_log10_generations"
stopifnot(nonconv %in% diag$model_id,
          round(diag$min_tail_ess[diag$model_id == nonconv]) == 385)

s4_tex_rows <- character(0)
s4_out <- list()
for (mod in cat_models$moderator) for (tv in c("log10_years", "log10_generations")) {
  id <- paste0(mod, "_plus_", tv)
  a  <- add[add$model_id == id, ]
  dg <- diag[diag$model_id == id, ]
  stopifnot(nrow(dg) == 1)
  levs <- levels(factor(dat_by_model[[cat_models$model_id[cat_models$moderator == mod]]][[mod]]))
  loc <- a[a$component == "location", ]; scl <- a[a$component == "scale", ]
  key <- function(x) ifelse(x$type == "slope", "slope", paste(x$level_a, x$level_b, sep = "|"))
  m <- merge(transform(loc, key = key(loc)), transform(scl, key = key(scl))[, c(
    "key", "estimate", "lower.CrI", "upper.CrI")], by = "key", suffixes = c("_loc", "_scl"))
  m$ord <- ifelse(m$type == "slope", 0,
                  match(m$level_b, levs) * 100 + match(m$level_a, levs))
  m <- m[order(m$ord), ]
  stopifnot(sum(m$type == "contrast") == choose(length(levs), 2))
  m$row_label <- ifelse(m$type == "slope",
                        paste0("Slope (", timevar_labels[[tv]], ")"),
                        paste0(mapply(function(x) if (is.na(x)) "" else lab(mod, x), m$level_a),
                               " $-$ ",
                               mapply(function(x) if (is.na(x)) "" else lab(mod, x), m$level_b)))
  flag <- if (id == nonconv) "$^{\\ddagger}$" else ""
  hdr <- paste0(moderator_labels[[mod]], " + ", timevar_labels[[tv]], flag)
  conv <- sprintf(paste0("$k=%s$, %d refs, %d spp; max $\\hat{R}$ %.3f; min bulk ESS %s; ",
                         "min tail ESS %s; divergences %d"),
                  fmt_int(a$k[1]), a$n_refs[1], a$n_species[1], dg$max_rhat,
                  fmt_int(floor(dg$min_bulk_ess)), fmt_int(floor(dg$min_tail_ess)),
                  dg$n_divergent)
  if (length(s4_tex_rows)) s4_tex_rows <- c(s4_tex_rows, "\\addlinespace")
  s4_tex_rows <- c(s4_tex_rows,
    paste0("\\multicolumn{3}{@{}l}{\\textbf{", hdr, "}} \\\\"),
    paste0("\\multicolumn{3}{@{}l}{\\quad\\textit{", conv, "}} \\\\"),
    paste0("\\quad ", m$row_label, " & ",
           fmt_ci(m$estimate_loc, m$lower.CrI_loc, m$upper.CrI_loc,
                  excl0(m$lower.CrI_loc, m$upper.CrI_loc)), " & ",
           fmt_ci(m$estimate_scl, m$lower.CrI_scl, m$upper.CrI_scl,
                  excl0(m$lower.CrI_scl, m$upper.CrI_scl)), " \\\\"))
  s4_out[[id]] <- data.frame(
    model_id = id, moderator = mod, timevar = tv, converged = id != nonconv,
    k = a$k[1], n_refs = a$n_refs[1], n_species = a$n_species[1],
    max_rhat = dg$max_rhat, min_bulk_ess = dg$min_bulk_ess,
    min_tail_ess = dg$min_tail_ess, n_divergent = dg$n_divergent,
    term = plain(m$row_label), type = m$type, minuend = m$level_a, subtrahend = m$level_b,
    loc_est = m$estimate_loc, loc_lo = m$lower.CrI_loc, loc_hi = m$upper.CrI_loc,
    scl_est = m$estimate_scl, scl_lo = m$lower.CrI_scl, scl_hi = m$upper.CrI_scl,
    stringsAsFactors = FALSE)
}
s4 <- do.call(rbind, s4_out); rownames(s4) <- NULL

s4_caption <- paste0(
  "\\textbf{Time-adjusted location--scale models (not preregistered).} ",
  "Each model added elapsed time ($\\log_{10}$ years or $\\log_{10}$ generations) ",
  "as an additive covariate to both the location and scale components of the ",
  "primary model for one categorical moderator. Values are posterior means with ",
  "95\\% equal-tailed credible intervals (CrI). Slopes are per one-unit increase ",
  "in $\\log_{10}$ elapsed time. Contrasts are posterior differences between ",
  "moderator coefficients, conditional on elapsed time; each contrast is the ",
  "first-named level minus the second-named level. Location values are on the ",
  "$\\ln M$ scale; scale values are on the scale of the natural logarithm of the ",
  "residual standard deviation, $\\ln\\sigma$. ",
  "The italic line under each model gives the number of contrasts ($k$), source ",
  "references (refs) and species (spp) analysed, and the convergence summary: ",
  "maximum $\\hat{R}$, minimum bulk and tail effective sample size (ESS), and the ",
  "number of divergent transitions. $^{\\dagger}$95\\% CrI excludes zero. ",
  "$^{\\ddagger}$This model did not meet the convergence criteria (minimum tail ",
  "ESS 385) and is shown for completeness only.")
s4_header <- "Term & Location [95\\% CrI] & Scale, $\\ln\\sigma$ [95\\% CrI] \\\\"
s4_tex <- c(ltx_begin("@{}lcc@{}", s4_caption, "tab:S4-time-adjusted", s4_header, 3),
            s4_tex_rows, ltx_end)

# ── 7. Cross-checks against values reported in the manuscript ─────────────────
xc <- tibble::tribble(
  ~model_id, ~component, ~minuend, ~subtrahend, ~est, ~lo, ~hi,
  "design_plus_log10_years",       "loc", "Synchronic", "Allochronic", 0.117, -0.058, 0.298,
  "design_plus_log10_generations", "loc", "Synchronic", "Allochronic", 0.111, -0.071, 0.290,
  "design_plus_log10_years",       "scl", "Synchronic", "Allochronic", -0.132, -0.195, -0.069,
  "design_plus_log10_generations", "scl", "Synchronic", "Allochronic", -0.056, -0.127, 0.015,
  "genphen_plus_log10_years",       "loc", "Phenotypic", "Genetic", -0.193, -0.336, -0.049,
  "genphen_plus_log10_generations", "loc", "Phenotypic", "Genetic", -0.183, -0.323, -0.048,
  "genphen_plus_log10_years",       "scl", "Phenotypic", "Genetic", 0.297, 0.239, 0.358,
  "disturbance_plus_log10_years",   "loc", "Landscape change", "Hunt_harv", 0.545, 0.195, 0.896,
  "disturbance_plus_log10_years",   "scl", "Other", "Introduction", 0.060, -0.013, 0.132
)
xc_res <- lapply(seq_len(nrow(xc)), function(i) {
  x <- xc[i, ]
  r <- s4[s4$model_id == x$model_id & s4$type == "contrast" &
            s4$minuend == x$minuend & s4$subtrahend == x$subtrahend, ]
  stopifnot(nrow(r) == 1)
  got <- unlist(r[paste0(x$component, c("_est", "_lo", "_hi"))])
  data.frame(model_id = x$model_id, component = x$component,
             contrast = paste(x$minuend, "-", x$subtrahend),
             manuscript = sprintf("%.3f [%.3f, %.3f]", x$est, x$lo, x$hi),
             table = sprintf("%.3f [%.3f, %.3f]", got[1], got[2], got[3]),
             match = isTRUE(all.equal(round(unname(got), 3), c(x$est, x$lo, x$hi))))
})
xc_res <- do.call(rbind, xc_res)
cat("\n=== Cross-check: Table S4 vs manuscript ===\n")
print(xc_res, row.names = FALSE)
if (!all(xc_res$match)) stop("Cross-check against manuscript values failed.")

# Primary-model checks: S3 vs manuscript text, S2 overall vs manuscript text.
p_chk <- data.frame(
  check = c("S3 design location Synchronic - Allochronic",
            "S3 design scale Synchronic - Allochronic",
            "S3 study type location Phenotypic - Genetic",
            "S2 overall location (2 dp)"),
  manuscript = c("0.250 [0.077, 0.418]", "-0.208 [-0.260, -0.157]",
                 "-0.210 [-0.355, -0.063]", "-0.55 [-1.41, 0.30]"),
  table = c(
    with(s3[s3$moderator == "design", ], sprintf("%.3f [%.3f, %.3f]", estimate_loc, lower_loc, upper_loc)),
    with(s3[s3$moderator == "design", ], sprintf("%.3f [%.3f, %.3f]", estimate_scl, lower_scl, upper_scl)),
    with(s3[s3$moderator == "genphen", ], sprintf("%.3f [%.3f, %.3f]", estimate_loc, lower_loc, upper_loc)),
    with(s2[s2$row == "Intercept" & s2$block == "Overall (intercept-only model)", ],
         sprintf("%.2f [%.2f, %.2f]", loc_est, loc_lo, loc_hi))))
p_chk$match <- p_chk$manuscript == p_chk$table
cat("\n=== Cross-check: Tables S2/S3 vs manuscript ===\n")
print(p_chk, row.names = FALSE)
if (!all(p_chk$match)) stop("Primary-model cross-check failed.")

# ── 8. Write outputs ──────────────────────────────────────────────────────────
s2_csv <- s2; s2_csv$row <- plain(s2_csv$row); s2_csv$block <- plain(s2_csv$block)
readr::write_csv(s2_csv, dir_out("tables", "esm_table_S2.csv"), na = "")
s3_csv <- s3[, c("block", "minuend", "subtrahend", "estimate_loc", "lower_loc",
                 "upper_loc", "estimate_scl", "lower_scl", "upper_scl",
                 "flipped_loc", "flipped_scl")]
s3_csv$block <- plain(s3_csv$block)
readr::write_csv(s3_csv, dir_out("tables", "esm_table_S3.csv"))
readr::write_csv(s4, dir_out("tables", "esm_table_S4.csv"), na = "")

check_balance(s2_tex, "S2"); check_balance(s3_tex, "S3"); check_balance(s4_tex, "S4")
write_tex(s2_tex, "esm_table_S2.tex")
write_tex(s3_tex, "esm_table_S3.tex")
write_tex(s4_tex, "esm_table_S4.tex")

cat("\n=== Row counts ===\n")
cat("S2:", sum(!s2$header), "data rows +", sum(s2$header), "block rows\n")
cat("S3:", nrow(s3), "contrast rows\n")
cat("S4:", nrow(s4), "data rows (", sum(s4$type == "slope"), "slope rows,",
    sum(s4$type == "contrast"), "contrast rows ) in", length(unique(s4$model_id)), "models\n")
cat("S3 rows re-oriented from the source CSV (location/scale):",
    sum(s3$flipped_loc), "/", sum(s3$flipped_scl), "\n")
