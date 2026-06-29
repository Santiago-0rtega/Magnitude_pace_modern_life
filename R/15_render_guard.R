# R/15_render_guard.R
# ─────────────────────────────────────────────────────────────────────────────
# Book render guard: the rendered book must READ cached .rds only — it must never
# recompute effect sizes, rebuild the database/phylogeny, or refit models.
#
# Sourcing this file turns on read-only mode for the *_or_read helpers
# (compute_or_read_lnm_safe, build_or_read_phylogeny, fit_or_read_model): if a
# required cache is absent they stop() with a clear message instead of silently
# launching hours of computation.
#
# This is OPT-IN per chapter and is NOT sourced by the model-fitting scripts on
# the compute server (run_m*.R, precompute_*.R), so those still fit/build as
# normal. An explicit recompute/refit = TRUE in a chapter also still overrides
# the guard, for the rare deliberate rebuild.
# ─────────────────────────────────────────────────────────────────────────────

options(pace.read_only = TRUE)

# Hard guard for a plain readRDS() call inside a chapter: confirm the cache exists
# before reading, otherwise stop with an actionable message.
pace_need <- function(path, what = "cache") {
  if (!file.exists(path)) {
    stop(sprintf(
      "[read-only render] Required %s not found:\n  %s\nThe book will not recompute. Build/sync this .rds first, then re-render.",
      what, path), call. = FALSE)
  }
  invisible(path)
}
