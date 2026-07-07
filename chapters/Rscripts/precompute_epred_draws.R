# precompute_epred_draws.R
# Run ONCE (re-run when a pending fit completes). Builds the small per-model
# draw caches in outputs/epred_draws/<id>.rds that back models/m*.qmd.
#
#   Rscript precompute_epred_draws.R          # build missing caches
#   FORCE=1 Rscript precompute_epred_draws.R  # rebuild all
#
# Heavy step (loads each ~400 MB fit + epreds); meant to run on totoro.

setwd(if (dir.exists("/home/ortegara/Documents/PACE"))
        "/home/ortegara/Documents/PACE" else getwd())
suppressMessages({ library(brms); library(tidybayes); library(dplyr); library(tibble) })
source(file.path("R", "14_epred_cache.R"))

force_all <- nzchar(Sys.getenv("FORCE"))
model_dir <- file.path("outputs", "models")
cache_dir <- file.path("outputs", "epred_draws")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# id, kind, candidate fit files (first existing wins) + matching moderator,
# label, predictors (interaction), draw settings.
specs <- list(
  list(id="m00", kind="intercept",   files="m00_ls_intercept_only",
       label="Intercept-only baseline", ndraws=NA),
  list(id="m01", kind="categorical",  files="m01_ls_disturbance",   mods="disturbance",
       label="Disturbance context", ndraws=1000),
  list(id="m02", kind="categorical",  files="m02_ls_design",        mods="design",
       label="Comparison design", ndraws=1000),
  list(id="m03", kind="continuous",   files="m03_ls_log10_years",   mods="log10_years",
       label="Elapsed time (log10 years)", ndraws=500, n_grid=100),
  list(id="m04", kind="continuous",   files="m04_ls_log10_generations", mods="log10_generations",
       label="Elapsed time (log10 generations)", ndraws=500, n_grid=100),
  list(id="m05", kind="categorical",  files="m05_ls_trait_type",    mods="trait_type",
       label="Trait type", ndraws=1000),
  # m06 (taxa) and m09 (data_type) are excluded from the primary results — see
  # R/06_model_registry.R (excluded_grid) and models/m06.qmd / m09.qmd for why.
  # No epred cache is built for them; those chapters don't call plot_cache().
  list(id="m07", kind="categorical",  files="m07_ls_genphen",       mods="genphen",
       label="Phenotypic vs genetic study", ndraws=1000),
  list(id="m08", kind="categorical",  files="m08_ls_env_change",    mods="env_change",
       label="Environmental-change context", ndraws=1000),
  list(id="m10", kind="categorical",  files="m10b_ls_transf_data_v2", mods="transf_data_v2",
       label="Transformation status", ndraws=1000),
  list(id="m11", kind="categorical",  files="m11_ls_data_scale",    mods="data_scale",
       label="Measurement scale", ndraws=1000),
  list(id="m12", kind="interaction",  files="m12_ls_years_x_disturbance",
       preds=c("log10_years","disturbance"),
       label="Elapsed time (log10 years)", ndraws=400, n_grid=80),
  list(id="m13", kind="interaction",  files="m13_ls_generations_x_disturbance",
       preds=c("log10_generations","disturbance"),
       label="Elapsed time (log10 generations)", ndraws=400, n_grid=80)
)

for (s in specs) {
  out_f <- file.path(cache_dir, paste0(s$id, ".rds"))
  if (file.exists(out_f) && !force_all) {
    message(s$id, ": cache exists — skip (FORCE=1 to rebuild)"); next
  }

  # pick first existing candidate fit
  hit <- which(file.exists(file.path(model_dir, paste0(s$files, ".rds"))))
  if (length(hit) == 0) {
    message(s$id, ": no fit yet (", paste(s$files, collapse=" / "), ") — skip"); next
  }
  i        <- hit[1]
  fit_file <- file.path(model_dir, paste0(s$files[i], ".rds"))
  mod      <- if (!is.null(s$mods)) s$mods[i] else NULL

  message("\n", s$id, ": building from ", basename(fit_file),
          if (!is.null(mod)) paste0("  [", mod, "]") else "")
  fit <- readRDS(fit_file)

  cache <- build_epred_cache(
    fit, id = s$id, kind = s$kind, label = s$label,
    moderator  = mod,
    predictors = s$preds,
    ndraws     = if (is.na(s$ndraws)) 1000 else s$ndraws,
    n_grid     = if (is.null(s$n_grid)) 100 else s$n_grid
  )
  cache$source_fit <- basename(fit_file)
  saveRDS(cache, out_f, compress = "xz")
  rm(fit, cache); gc(verbose = FALSE)
  message("  saved ", out_f, "  (",
          round(file.size(out_f) / 1024^2, 1), " MB)")
}

message("\nAll done.")
