# Build a brms bf() for a single location-scale model.
# moderator: bare string name of the predictor column.
# has_phylogeny: whether to include the gr(sp_ncbi, cov = A) term.
# has_sys_id: whether to include (1 | sys_id) for sensitivity analyses.
build_ls_formula <- function(moderator,
                             has_phylogeny = TRUE,
                             has_sys_id    = FALSE) {

  loc_terms <- paste0(
    "yi_lnM_safe ~ ", moderator,
    " + (1 | ref_id)",
    if (has_sys_id)    " + (1 | sys_id)"              else "",
    if (has_phylogeny) " + (1 | gr(sp_ncbi, cov = A))" else "",
    " + (1 | gr(es_id_model, cov = V))"
  )

  scl_terms <- paste0("sigma ~ ", moderator)

  brms::bf(
    stats::as.formula(loc_terms),
    stats::as.formula(scl_terms)
  )
}

# ── Prior helpers ─────────────────────────────────────────────────────────────

# Function 1: default priors for all parameters.
# Returns the raw brmsprior table from default_prior().
get_default_priors <- function(formula, dat_model, V, A = NULL) {
  data2 <- list(V = V)
  if (!is.null(A)) data2$A <- A

  brms::default_prior(
    formula,
    data   = dat_model,
    data2  = data2,
    family = gaussian()
  )
}

# Function 2: set constant(1) for the es_id_model SD.
#
# Strategy: remove ALL existing SD rows for the es_id_model group from the
# prior table, then replace them with a single fresh prior_string() row that
# targets the coefficient level (coef = "Intercept") directly.
#
# Why not just overwrite in place?
# default_prior() returns two rows for each random-effects group SD:
#   coef = ""          — group-level ("global for group"), less specific
#   coef = "Intercept" — coefficient-level, MORE SPECIFIC, takes precedence
#
# Modifying only the group-level row causes brms to warn that its
# more-specific Intercept-level prior (still at the default) overrides ours.
# Overwriting both in place is fragile because rbind on brmsprior objects
# can drop the class. The safe approach is to strip and rebuild.
set_esid_prior_constant <- function(pr) {

  # Drop every SD row belonging to the es_id_model group
  pr_clean <- pr[!(pr$class == "sd" & pr$group == "es_id_model"), ]

  # Build a fresh constant(1) prior targeted at the Intercept (coefficient) level.
  # prior_string() produces a properly-classed brmsprior object.
  pr_const <- brms::prior_string(
    "constant(1)",
    class = "sd",
    group = "es_id_model",
    coef  = "Intercept"
  )

  # Combine: pr_clean retains the brmsprior class, so rbind preserves it.
  rbind(pr_clean, pr_const)
}

# Function 3: verify that constant(1) is correctly set before fitting.
# Prints the es_id_model SD rows and returns TRUE/FALSE invisibly.
verify_esid_prior <- function(pr, model_name = "") {
  esid_rows <- pr[pr$class == "sd" & pr$group == "es_id_model", ]

  prefix <- if (nzchar(model_name)) paste0("[", model_name, "] ") else ""

  if (nrow(esid_rows) == 0) {
    message(prefix, "PRIOR CHECK FAIL — no SD row found for group 'es_id_model'.")
    return(invisible(FALSE))
  }

  ok <- any(esid_rows$coef == "Intercept" & esid_rows$prior == "constant(1)")

  status <- if (ok) "OK" else "WARN — constant(1) NOT at coef=Intercept level"

  message(
    prefix, "Prior check: ", status, "\n",
    "  es_id_model SD rows in prior table:\n",
    paste0(
      "    coef='", esid_rows$coef, "'  prior='", esid_rows$prior, "'",
      collapse = "\n"
    )
  )

  invisible(ok)
}

# Convenience wrapper: build full priors with constant(1) already applied.
build_ls_priors <- function(formula, dat_model, V, A = NULL) {
  pr <- get_default_priors(formula, dat_model, V, A)
  pr <- set_esid_prior_constant(pr)
  pr
}
