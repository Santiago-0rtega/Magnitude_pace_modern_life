# ─────────────────────────────────────────────────────────────────────────────
# run_sensitivity.R  —  Protocol §3.7 sensitivity analyses
#
# Fits the three coded sensitivity variants for every included moderator model
# (m01–m05, m07, m08 — i.e. book result-chapters except the m12 interaction and
# the excluded m06/m09):
#
#   Sensitivity 1  — minimum total sample size (n1 + n2 >= 40)      [protocol §3.7 #1]
#   Sensitivity 2  — no-phylogeny fallback (drop the sp_ncbi/A term)
#   Sensitivity 3  — add (1 | sys_id) random effect                [protocol §3.7 #4]
#
# 7 moderators × 3 variants = 21 models. Models cached to
# outputs/models/sensitivity/ ; summary CSVs to outputs/tables/sensitivity/.
# Re-running skips any model whose .rds already exists (refit_* = FALSE).
# ─────────────────────────────────────────────────────────────────────────────

setwd(here::here())
library(here)

source(here::here("chapters", "Rscripts", "00_packages_fit.R"))
source(here::here("chapters", "Rscripts", "01_paths.R"))
source(here::here("chapters", "Rscripts", "05_phylogeny.R"))
source(here::here("chapters", "Rscripts", "06_model_registry.R"))
source(here::here("chapters", "Rscripts", "07_model_formulas.R"))
source(here::here("chapters", "Rscripts", "08_fit_or_read_model.R"))
source(here::here("chapters", "Rscripts", "09_model_summaries.R"))
source(here::here("chapters", "Rscripts", "10_model_diagnostics.R"))

# Never silently skip because of the read-only render guard.
options(pace.read_only = FALSE)

# Fit if missing; set TRUE to force refit of already-cached sensitivity models.
refit_sensitivity <- FALSE

# Moderators = included book models except m12 (interaction) and m06/m09 (excluded).
#   disturbance=m01, design=m02, log10_years=m03, log10_generations=m04,
#   trait_type=m05, genphen=m07, env_change=m08
priority_mods <- c("disturbance", "design", "log10_years",
                   "log10_generations", "trait_type", "genphen", "env_change")

model_dir <- dir_out("models", "sensitivity")
table_dir <- dir_out("tables", "sensitivity")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)

# ── Data ─────────────────────────────────────────────────────────────────────
dat_es <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
stamp("Total contrasts: ", nrow(dat_es))

A_full   <- tryCatch(readRDS(dir_out("phylogeny", "proceed_A_matrix.rds")),
                     error = function(e) NULL)
name_map <- tryCatch(readRDS(dir_out("phylogeny", "proceed_name_map.rds")),
                     error = function(e) NULL)

if (!is.null(name_map)) {
  dat_es <- apply_phylo_name_map(dat_es, name_map)
} else {
  dat_es$sp_ncbi_canonical <- dat_es$sp_ncbi
}

# ═════════════════════════════════════════════════════════════════════════════
# Sensitivity 1: minimum total sample size (n1 + n2 >= 40)
# ═════════════════════════════════════════════════════════════════════════════
stamp("===== Sensitivity 1: N >= 40 =====")
dat_n40 <- dat_es |> dplyr::filter(n_total >= 40)
stamp("Contrasts with N >= 40: ", nrow(dat_n40))

sens1_results <- list()
for (moderator in priority_mods) {
  if (!moderator %in% names(dat_n40)) next
  dat_mod <- dat_n40[!is.na(dat_n40[[moderator]]), ]
  if (nrow(dat_mod) < 10) { stamp("skip sens1 ", moderator, ": too few rows"); next }

  dat_mod <- dat_mod |> dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
  if (is.factor(dat_mod[[moderator]])) dat_mod[[moderator]] <- droplevels(dat_mod[[moderator]])

  phylo     <- prepare_phylo_and_data(dat_mod, A_full, label = paste0("sens1_n40_", moderator))
  dat_mod   <- phylo$dat_model
  A_mod     <- phylo$A_mod
  V         <- phylo$V
  has_phylo <- phylo$has_phylo

  formula <- build_ls_formula(moderator, has_phylogeny = has_phylo)
  priors  <- build_ls_priors(formula, dat_mod, V, A = A_mod)
  mname   <- paste0("sens1_n40_", moderator)
  if (!verify_esid_prior(priors, model_name = mname)) { stamp("skip ", mname, ": prior check failed"); next }

  stamp("fitting ", mname, " (rows=", nrow(dat_mod), ", phylo=", has_phylo, ")")
  fit <- tryCatch(
    fit_or_read_model(mname,
                      fit_fun   = function() fit_ls_model(dat_mod, formula, priors, V = V, A = A_mod,
                                                          mcmc_args = default_mcmc_args),
                      model_dir = model_dir, refit = refit_sensitivity),
    error = function(e) { stamp("FAILED ", mname, ": ", conditionMessage(e)); NULL })

  if (!is.null(fit)) sens1_results[[moderator]] <- list(
    fe   = extract_fixed_effects(fit, mname, moderator),
    diag = tryCatch(extract_diagnostics(fit, mname, moderator), error = function(e) NULL))
}
if (length(sens1_results) > 0) {
  readr::write_csv(purrr::map_dfr(sens1_results, ~ dplyr::filter(.x$fe, submodel == "location")),
                   file.path(table_dir, "sens1_n40_location.csv"))
  readr::write_csv(purrr::map_dfr(sens1_results, ~ dplyr::filter(.x$fe, submodel == "scale")),
                   file.path(table_dir, "sens1_n40_scale.csv"))
  readr::write_csv(purrr::map_dfr(purrr::compact(purrr::map(sens1_results, "diag")), ~ .x),
                   file.path(table_dir, "sens1_n40_diagnostics.csv"))
}

# ═════════════════════════════════════════════════════════════════════════════
# Sensitivity 2: no-phylogeny fallback
# ═════════════════════════════════════════════════════════════════════════════
stamp("===== Sensitivity 2: no-phylogeny =====")
sens2_results <- list()
for (moderator in priority_mods) {
  if (!moderator %in% names(dat_es)) next
  dat_mod <- dat_es[!is.na(dat_es[[moderator]]), ]
  if (nrow(dat_mod) < 10) next

  dat_mod <- dat_mod |> dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
  if (is.factor(dat_mod[[moderator]])) dat_mod[[moderator]] <- droplevels(dat_mod[[moderator]])

  V <- diag(dat_mod$vi_lnM_safe)
  rownames(V) <- colnames(V) <- as.character(dat_mod$es_id_model)

  formula <- build_ls_formula(moderator, has_phylogeny = FALSE)
  priors  <- build_ls_priors(formula, dat_mod, V, A = NULL)
  mname   <- paste0("sens2_nophylo_", moderator)
  if (!verify_esid_prior(priors, model_name = mname)) { stamp("skip ", mname, ": prior check failed"); next }

  stamp("fitting ", mname, " (rows=", nrow(dat_mod), ")")
  fit <- tryCatch(
    fit_or_read_model(mname,
                      fit_fun   = function() fit_ls_model(dat_mod, formula, priors, V = V, A = NULL,
                                                          mcmc_args = default_mcmc_args),
                      model_dir = model_dir, refit = refit_sensitivity),
    error = function(e) { stamp("FAILED ", mname, ": ", conditionMessage(e)); NULL })

  if (!is.null(fit)) sens2_results[[moderator]] <- list(
    fe   = extract_fixed_effects(fit, mname, moderator),
    diag = tryCatch(extract_diagnostics(fit, mname, moderator), error = function(e) NULL))
}
if (length(sens2_results) > 0) {
  readr::write_csv(purrr::map_dfr(sens2_results, ~ dplyr::filter(.x$fe, submodel == "location")),
                   file.path(table_dir, "sens2_nophylo_location.csv"))
  readr::write_csv(purrr::map_dfr(sens2_results, ~ dplyr::filter(.x$fe, submodel == "scale")),
                   file.path(table_dir, "sens2_nophylo_scale.csv"))
  readr::write_csv(purrr::map_dfr(purrr::compact(purrr::map(sens2_results, "diag")), ~ .x),
                   file.path(table_dir, "sens2_nophylo_diagnostics.csv"))
}

# ═════════════════════════════════════════════════════════════════════════════
# Sensitivity 3: sys_id random effect
# ═════════════════════════════════════════════════════════════════════════════
stamp("===== Sensitivity 3: sys_id =====")
sens3_results <- list()
for (moderator in priority_mods) {
  if (!moderator %in% names(dat_es)) next
  dat_mod <- dat_es[!is.na(dat_es[[moderator]]) & !is.na(dat_es$sys_id), ]
  if (nrow(dat_mod) < 10) next

  dat_mod <- dat_mod |>
    dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())),
                  sys_id      = droplevels(factor(sys_id)))
  if (is.factor(dat_mod[[moderator]])) dat_mod[[moderator]] <- droplevels(dat_mod[[moderator]])

  phylo     <- prepare_phylo_and_data(dat_mod, A_full, label = paste0("sens3_sysid_", moderator))
  dat_mod   <- phylo$dat_model
  A_mod     <- phylo$A_mod
  V         <- phylo$V
  has_phylo <- phylo$has_phylo

  formula <- build_ls_formula(moderator, has_phylogeny = has_phylo, has_sys_id = TRUE)
  priors  <- build_ls_priors(formula, dat_mod, V, A = A_mod)
  mname   <- paste0("sens3_sysid_", moderator)
  if (!verify_esid_prior(priors, model_name = mname)) { stamp("skip ", mname, ": prior check failed"); next }

  stamp("fitting ", mname, " (rows=", nrow(dat_mod), ", phylo=", has_phylo, ")")
  fit <- tryCatch(
    fit_or_read_model(mname,
                      fit_fun   = function() fit_ls_model(dat_mod, formula, priors, V = V, A = A_mod,
                                                          mcmc_args = default_mcmc_args),
                      model_dir = model_dir, refit = refit_sensitivity),
    error = function(e) { stamp("FAILED ", mname, ": ", conditionMessage(e)); NULL })

  if (!is.null(fit)) {
    diag_row <- tryCatch(extract_diagnostics(fit, mname, moderator), error = function(e) NULL)
    if (!is.null(diag_row) && !is.na(diag_row$n_divergent) && diag_row$n_divergent > 20)
      stamp("WARNING ", mname, ": ", diag_row$n_divergent,
            " divergent transitions — may be weakly identifiable with sys_id.")
    sens3_results[[moderator]] <- list(
      fe   = extract_fixed_effects(fit, mname, moderator),
      diag = diag_row)
  }
}
if (length(sens3_results) > 0) {
  readr::write_csv(purrr::map_dfr(sens3_results, ~ dplyr::filter(.x$fe, submodel == "location")),
                   file.path(table_dir, "sens3_sysid_location.csv"))
  readr::write_csv(purrr::map_dfr(sens3_results, ~ dplyr::filter(.x$fe, submodel == "scale")),
                   file.path(table_dir, "sens3_sysid_scale.csv"))
  readr::write_csv(purrr::map_dfr(purrr::compact(purrr::map(sens3_results, "diag")), ~ .x),
                   file.path(table_dir, "sens3_sysid_diagnostics.csv"))
}

stamp("===== DONE =====")
stamp("models  -> ", model_dir)
stamp("tables  -> ", table_dir)
print(list.files(model_dir, pattern = "\\.rds$"))
