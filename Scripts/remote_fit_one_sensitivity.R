# ─────────────────────────────────────────────────────────────────────────────
# fit_one_sensitivity.R  <variant> <moderator>
#
# Fits exactly ONE sensitivity model so the 21-model grid can be run as
# independent parallel jobs (capped concurrency via xargs -P).
#
#   variant   ∈ { sens1_n40 , sens2_nophylo , sens3_sysid }
#   moderator ∈ { disturbance, design, log10_years, log10_generations,
#                 trait_type, genphen, env_change }
#
# Saves  outputs/models/sensitivity/<variant>_<moderator>.rds  and per-model
# parts  outputs/tables/sensitivity/parts/<variant>_<moderator>_{fe,diag}.csv
# Skips fitting if the .rds already exists (idempotent / resumable).
# ─────────────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: fit_one_sensitivity.R <variant> <moderator>")
variant   <- args[[1]]
moderator <- args[[2]]
mname     <- paste0(variant, "_", moderator)

setwd(here::here())
suppressMessages({
  library(here)
  for (f in c("00_packages_fit","01_paths","05_phylogeny","06_model_registry",
              "07_model_formulas","08_fit_or_read_model","09_model_summaries",
              "10_model_diagnostics"))
    source(here::here("Scripts", paste0(f, ".R")))
})
options(pace.read_only = FALSE)
refit_sensitivity <- identical(tolower(Sys.getenv("REFIT_SENSITIVITY", "false")), "true")

model_dir <- dir_out("models", "sensitivity")
part_dir  <- dir_out("tables", "sensitivity", "parts")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(part_dir,  recursive = TRUE, showWarnings = FALSE)

stamp <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  [", mname, "] ", ...)

# ── Data ─────────────────────────────────────────────────────────────────────
dat_es   <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
A_full   <- tryCatch(readRDS(dir_out("phylogeny", "proceed_A_matrix.rds")),
                     error = function(e) NULL)
name_map <- tryCatch(readRDS(dir_out("phylogeny", "proceed_name_map.rds")),
                     error = function(e) NULL)
if (!is.null(name_map)) {
  dat_es <- apply_phylo_name_map(dat_es, name_map)
} else {
  dat_es$sp_ncbi_canonical <- dat_es$sp_ncbi
}

if (!moderator %in% names(dat_es)) stop("moderator not in data: ", moderator)

# ── Build the one (data, formula, priors) triple for this variant ────────────
if (variant == "sens1_n40") {
  dat_mod <- dat_es |> dplyr::filter(n_total >= 40)
  dat_mod <- dat_mod[!is.na(dat_mod[[moderator]]), ]
  stopifnot(nrow(dat_mod) >= 10)
  dat_mod <- dat_mod |> dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
  if (is.factor(dat_mod[[moderator]])) dat_mod[[moderator]] <- droplevels(dat_mod[[moderator]])
  phylo   <- prepare_phylo_and_data(dat_mod, A_full, label = mname)
  dat_mod <- phylo$dat_model; A_mod <- phylo$A_mod; V <- phylo$V; has_phylo <- phylo$has_phylo
  formula <- build_ls_formula(moderator, has_phylogeny = has_phylo)

} else if (variant == "sens2_nophylo") {
  dat_mod <- dat_es[!is.na(dat_es[[moderator]]), ]
  stopifnot(nrow(dat_mod) >= 10)
  dat_mod <- dat_mod |> dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
  if (is.factor(dat_mod[[moderator]])) dat_mod[[moderator]] <- droplevels(dat_mod[[moderator]])
  V <- diag(dat_mod$vi_lnM_safe); rownames(V) <- colnames(V) <- as.character(dat_mod$es_id_model)
  A_mod <- NULL; has_phylo <- FALSE
  formula <- build_ls_formula(moderator, has_phylogeny = FALSE)

} else if (variant == "sens3_sysid") {
  dat_mod <- dat_es[!is.na(dat_es[[moderator]]) & !is.na(dat_es$sys_id), ]
  stopifnot(nrow(dat_mod) >= 10)
  dat_mod <- dat_mod |>
    dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())),
                  sys_id      = droplevels(factor(sys_id)))
  if (is.factor(dat_mod[[moderator]])) dat_mod[[moderator]] <- droplevels(dat_mod[[moderator]])
  phylo   <- prepare_phylo_and_data(dat_mod, A_full, label = mname)
  dat_mod <- phylo$dat_model; A_mod <- phylo$A_mod; V <- phylo$V; has_phylo <- phylo$has_phylo
  formula <- build_ls_formula(moderator, has_phylogeny = has_phylo, has_sys_id = TRUE)

} else stop("unknown variant: ", variant)

priors <- build_ls_priors(formula, dat_mod, V, A = A_mod)
if (!verify_esid_prior(priors, model_name = mname)) stop("prior check failed: ", mname)

stamp("fitting (rows=", nrow(dat_mod), ", phylo=", has_phylo, ")")
t0  <- Sys.time()
fit <- fit_or_read_model(
  mname,
  fit_fun   = function() fit_ls_model(dat_mod, formula, priors, V = V, A = A_mod,
                                      mcmc_args = default_mcmc_args),
  model_dir = model_dir, refit = refit_sensitivity)
stamp("done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

# ── Per-model summary parts ──────────────────────────────────────────────────
fe <- extract_fixed_effects(fit, mname, moderator)
readr::write_csv(fe, file.path(part_dir, paste0(mname, "_fe.csv")))
diag <- tryCatch(extract_diagnostics(fit, mname, moderator), error = function(e) NULL)
if (!is.null(diag)) readr::write_csv(diag, file.path(part_dir, paste0(mname, "_diag.csv")))
if (!is.null(diag) && !is.na(diag$n_divergent) && diag$n_divergent > 20)
  stamp("WARNING: ", diag$n_divergent, " divergent transitions.")
stamp("OK")
