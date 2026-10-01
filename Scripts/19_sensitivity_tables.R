# ─────────────────────────────────────────────────────────────────────────────
# 19_sensitivity_tables.R
# Read-side helpers for the sensitivity appendix chapters. They consume the
# small CSVs synced from totoro:
#   Rdata/tables/sensitivity/<variant>_{location,scale,diagnostics}.csv
#   Rdata/tables/primary/primary_fixef.csv           (matched primary fits)
# No brms fit is ever loaded at render.
#
# Every sensitivity model is a raw brms location–scale fit; the *_location /
# *_scale CSVs hold its fixed effects on the model's default (treatment) coding:
#   - Intercept        = reference-level mean (location) / log-SD (scale)
#   - other rows       = difference from the reference level
# This differs from the main model chapters, which report emmeans marginal
# means; here we compare raw coefficients against the matched primary fit so the
# two columns are strictly like-for-like.
# ─────────────────────────────────────────────────────────────────────────────

suppressMessages({
  library(dplyr)
  library(readr)
  library(kableExtra)
})

# Variant metadata -------------------------------------------------------------
sens_variant_meta <- function(variant) {
  switch(variant,
    sens1_n40     = list(short = "N ≥ 40",
                         title = "Minimum total sample size (n₁ + n₂ ≥ 40)"),
    sens2_nophylo = list(short = "no phylogeny",
                         title = "No phylogenetic random effect"),
    sens3_sysid   = list(short = "sys_id",
                         title = "Study-system (sys_id) random effect"),
    stop("unknown variant: ", variant))
}

# Moderator display labels + canonical order -----------------------------------
.sens_mod_labels <- c(
  disturbance       = "Disturbance context (m01)",
  design            = "Comparison design (m02)",
  log10_years       = "Elapsed time — log₁₀ years (m03)",
  log10_generations = "Elapsed time — log₁₀ generations (m04)",
  trait_type        = "Trait type (m05)",
  genphen           = "Phenotypic vs genetic study (m07)"
)
sens_moderators <- names(.sens_mod_labels)

# Moderators excluded from reporting. m08 (env_change) failed convergence in
# the primary grid (max Rhat 1.0143) and is excluded from the primary results,
# so its sensitivity variants are excluded too. The fitted CSV rows are kept on
# disk but dropped at read time.
.sens_excluded_moderators <- c("env_change")

# Readers ----------------------------------------------------------------------
.sens_read <- function(variant, kind) {
  f <- here::here("Rdata", "tables", "sensitivity", paste0(variant, "_", kind, ".csv"))
  if (!file.exists(f)) return(NULL)
  readr::read_csv(f, show_col_types = FALSE) |>
    dplyr::filter(!moderator %in% .sens_excluded_moderators)
}
.primary_read <- function() {
  f <- here::here("Rdata", "tables", "primary", "primary_fixef.csv")
  if (!file.exists(f)) return(NULL)
  readr::read_csv(f, show_col_types = FALSE) |>
    dplyr::filter(!moderator %in% .sens_excluded_moderators)
}

# Formatting helpers -----------------------------------------------------------
.excl0   <- function(lo, hi) !is.na(lo) & !is.na(hi) & ((lo > 0 & hi > 0) | (lo < 0 & hi < 0))
.fmt_ci  <- function(est, lo, hi) ifelse(is.na(est), "—",
                                         sprintf("%.3f [%.3f, %.3f]", est, lo, hi))

.pretty_term <- function(term, moderator) {
  base <- sub("^sigma_", "", term)
  vapply(seq_along(base), function(i) {
    b <- base[i]
    if (b == "Intercept") return("Intercept (ref.)")
    if (b == moderator)   return("slope")
    lev <- sub(paste0("^", moderator), "", b)
    if (lev == "") lev <- b
    lev
  }, character(1))
}

# Per-moderator comparison table: sensitivity vs matched primary fit -----------
sens_moderator_table <- function(variant, moderator) {
  meta <- sens_variant_meta(variant)
  loc  <- .sens_read(variant, "location")
  scl  <- .sens_read(variant, "scale")
  sens <- dplyr::bind_rows(loc, scl) |> dplyr::filter(moderator == !!moderator)
  if (is.null(sens) || nrow(sens) == 0)
    return(kableExtra::kbl(data.frame(Note = "No sensitivity result found."),
                           caption = .sens_mod_labels[[moderator]]))

  prim <- .primary_read()
  prim <- if (is.null(prim)) NULL else
    prim |> dplyr::filter(moderator == !!moderator) |>
      dplyr::select(submodel, term, p_est = estimate, p_lo = q2_5, p_hi = q97_5)

  tab <- sens |>
    dplyr::transmute(
      submodel, term,
      s_est = estimate, s_lo = q2_5, s_hi = q97_5)
  if (!is.null(prim)) tab <- dplyr::left_join(tab, prim, by = c("submodel", "term"))
  else tab <- tab |> dplyr::mutate(p_est = NA_real_, p_lo = NA_real_, p_hi = NA_real_)

  tab <- tab |>
    dplyr::mutate(
      s_excl = .excl0(s_lo, s_hi),
      p_excl = .excl0(p_lo, p_hi),
      Stability = dplyr::case_when(
        is.na(p_est)                               ~ "level absent in primary",
        s_excl == p_excl & sign(s_est) == sign(p_est) ~ "consistent",
        TRUE                                       ~ "DIFFERS"),
      Sub          = ifelse(submodel == "location", "Location (mean lnM)", "Scale (log-SD)"),
      Term         = .pretty_term(term, moderator),
      `Sensitivity β [95% CrI]` = .fmt_ci(s_est, s_lo, s_hi),
      `Primary β [95% CrI]`     = .fmt_ci(p_est, p_lo, p_hi)) |>
    dplyr::arrange(dplyr::desc(submodel == "location"))  # location block first

  disp <- tab |>
    dplyr::select(`Sub-model` = Sub, Term,
                  `Sensitivity β [95% CrI]`, `Primary β [95% CrI]`, Stability)

  kb <- kableExtra::kbl(disp, caption = .sens_mod_labels[[moderator]], booktabs = TRUE) |>
    kableExtra::kable_styling(full_width = FALSE,
                              bootstrap_options = c("striped", "condensed"))
  kb
}

# Convergence-diagnostics table for one variant --------------------------------
sens_diagnostics_table <- function(variant) {
  d <- .sens_read(variant, "diagnostics")
  if (is.null(d) || nrow(d) == 0)
    return(kableExtra::kbl(data.frame(Note = "No diagnostics found.")))
  d <- d |>
    dplyr::mutate(Moderator = .sens_mod_labels[as.character(moderator)]) |>
    dplyr::transmute(
      Moderator,
      `Max R̂`      = round(max_rhat, 4),
      `Min bulk ESS`     = round(min_bulk_ess),
      `Min tail ESS`     = round(min_tail_ess),
      Divergences        = n_divergent,
      `Max-treedepth hits` = max_treedepth_hits) |>
    dplyr::arrange(dplyr::desc(Divergences))
  kb <- kableExtra::kbl(d, caption = "Convergence diagnostics (4 chains, 2000 post-warmup draws each).",
                        booktabs = TRUE) |>
    kableExtra::kable_styling(full_width = FALSE,
                              bootstrap_options = c("striped", "condensed"))
  kb
}

# Render-time diagnostic note, generated from the latest synced table so prose
# cannot retain divergence/Rhat values from an older rerun.
sens_diagnostic_flags <- function(variant) {
  d <- .sens_read(variant, "diagnostics")
  if (is.null(d) || nrow(d) == 0) return("No diagnostic results are available.")
  diagnostic_labels <- c(
    disturbance = "Disturbance context (m01)",
    design = "Comparison design (m02)",
    log10_years = "Elapsed time in log10 years (m03)",
    log10_generations = "Elapsed time in log10 generations (m04)",
    trait_type = "Trait type (m05)",
    genphen = "Phenotypic vs genetic study (m07)"
  )
  # Manuscript convergence criteria: max Rhat <= 1.01, bulk and tail ESS >= 400,
  # and no divergent transitions.
  flagged <- d |>
    dplyr::filter(n_divergent > 0 | max_rhat > 1.01 |
                    min_bulk_ess < 400 | min_tail_ess < 400) |>
    dplyr::mutate(
      label = unname(diagnostic_labels[as.character(moderator)]),
      detail = sprintf("%s: max Rhat = %.4f, min bulk/tail ESS = %.0f/%.0f, divergences = %d",
                       label, max_rhat, min_bulk_ess, min_tail_ess, n_divergent)
    )
  if (nrow(flagged) == 0) {
    return(sprintf("All %d models had max Rhat <= 1.01, bulk and tail ESS >= 400, and zero divergences (largest max Rhat = %.4f).",
                   nrow(d), max(d$max_rhat, na.rm = TRUE)))
  }
  paste0(
    "Models requiring convergence caution: ",
    paste(flagged$detail, collapse = "; "),
    ". See the table above for ESS and treedepth diagnostics."
  )
}

# Emit every moderator's comparison table (for a results='asis' chunk) ---------
sens_all_moderator_tables <- function(variant) {
  out <- lapply(sens_moderators, function(m) {
    knitr::knit_child(
      text = c(sprintf("### %s\n", .sens_mod_labels[[m]]),
               "```{r}",
               sprintf("sens_moderator_table('%s', '%s')", variant, m),
               "```\n"),
      quiet = TRUE, envir = environment())
  })
  cat(unlist(out), sep = "\n")
}

# One-line count of models whose conclusion differs from primary ---------------
sens_stability_summary <- function(variant) {
  loc <- .sens_read(variant, "location"); scl <- .sens_read(variant, "scale")
  prim <- .primary_read()
  if (is.null(prim)) return("Primary comparison unavailable.")
  s <- dplyr::bind_rows(loc, scl) |>
    dplyr::select(moderator, submodel, term, s_est = estimate, s_lo = q2_5, s_hi = q97_5)
  p <- prim |> dplyr::select(moderator, submodel, term, p_lo = q2_5, p_hi = q97_5, p_est = estimate)
  j <- dplyr::left_join(s, p, by = c("moderator", "submodel", "term")) |>
    dplyr::filter(term != "Intercept", term != "sigma_Intercept", !is.na(p_est)) |>
    dplyr::mutate(differs = .excl0(s_lo, s_hi) != .excl0(p_lo, p_hi) | sign(s_est) != sign(p_est))
  sprintf("%d of %d moderator coefficients change their credible-interval conclusion relative to the primary fit.",
          sum(j$differs, na.rm = TRUE), nrow(j))
}

# Book-level synthesis across all sensitivity variants ------------------------
sens_variant_summary_data <- function() {
  variants <- c("sens1_n40", "sens2_nophylo", "sens3_sysid")

  purrr::map_dfr(variants, function(variant) {
    loc  <- .sens_read(variant, "location")
    scl  <- .sens_read(variant, "scale")
    diag <- .sens_read(variant, "diagnostics")
    prim <- .primary_read()

    if (is.null(loc) || is.null(scl) || is.null(diag)) {
      return(tibble::tibble(
        Analysis = sens_variant_meta(variant)$title,
        `Models available` = if (is.null(diag)) 0L else nrow(diag),
        `Changed coefficients` = NA_character_,
        `Max Rhat` = NA_real_, Divergences = NA_integer_, Status = "incomplete"
      ))
    }

    sens <- dplyr::bind_rows(loc, scl) |>
      dplyr::select(moderator, submodel, term,
                    s_est = estimate, s_lo = q2_5, s_hi = q97_5)
    joined <- sens |>
      dplyr::left_join(
        prim |>
          dplyr::select(moderator, submodel, term,
                        p_est = estimate, p_lo = q2_5, p_hi = q97_5),
        by = c("moderator", "submodel", "term")
      ) |>
      dplyr::filter(term != "Intercept", term != "sigma_Intercept", !is.na(p_est)) |>
      dplyr::mutate(
        differs = .excl0(s_lo, s_hi) != .excl0(p_lo, p_hi) |
          sign(s_est) != sign(p_est)
      )

    n_changed <- sum(joined$differs, na.rm = TRUE)
    n_compared <- nrow(joined)
    n_divergent <- sum(diag$n_divergent, na.rm = TRUE)

    tibble::tibble(
      Analysis = sens_variant_meta(variant)$title,
      `Models available` = nrow(diag),
      `Changed coefficients` = sprintf("%d of %d", n_changed, n_compared),
      `Max Rhat` = max(diag$max_rhat, na.rm = TRUE),
      Divergences = n_divergent,
      Status = if (nrow(diag) == length(sens_moderators)) "complete" else "incomplete"
    )
  })
}

sens_variant_summary_table <- function() {
  out <- sens_variant_summary_data() |>
    dplyr::mutate(`Max Rhat` = round(`Max Rhat`, 3))

  kableExtra::kbl(
    out,
    caption = "Sensitivity-model availability, stability, and diagnostics.",
    booktabs = TRUE
  ) |>
    kableExtra::kable_styling(
      full_width = FALSE,
      bootstrap_options = c("striped", "condensed")
    )
}
