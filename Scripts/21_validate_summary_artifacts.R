# Reject legacy numerical artifacts before Quarto executes any chapter.
# Specification v2 is posterior means with equal-tail 95% credible intervals.

spec_version <- 2L
need <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  readRDS(path)
}

summary_path <- file.path("Rdata", "summaries", "emmeans_contrasts_cache.rds")
summaries <- need(summary_path, "moderator summary cache")
if (!identical(attr(summaries, "summary_spec")$version, spec_version)) {
  stop("Legacy moderator summary cache rejected. Run precompute_emmeans_contrasts.R and sync the v2 cache.",
       call. = FALSE)
}

epred_ids <- c("m00", "m01", "m02", "m03", "m04", "m05", "m07", "m08", "m11")
for (id in epred_ids) {
  cache <- need(file.path("Rdata", "epred_draws", paste0(id, ".rds")),
                paste(id, "prediction cache"))
  if (!identical(cache$summary_spec$version, spec_version)) {
    stop("Legacy prediction cache rejected for ", id,
         ". Run precompute_epred_draws.R with FORCE=1 and sync all v2 caches.",
         call. = FALSE)
  }
}

manifest <- need(
  file.path("Rdata", "figures", "publication", "orchard", "orchard_manifest.rds"),
  "orchard figure manifest"
)
if (!identical(manifest$summary_spec$version, spec_version)) {
  stop("Legacy orchard figures rejected. Regenerate them from the v2 caches.",
       call. = FALSE)
}
if (!setequal(manifest$model_ids, epred_ids)) {
  stop("Orchard manifest does not cover every active moderator model.", call. = FALSE)
}
current_summary_md5 <- unname(tools::md5sum(summary_path))
if (!identical(manifest$summary_cache_md5, current_summary_md5)) {
  stop("Orchard figures are older than the active moderator summary cache.", call. = FALSE)
}
epred_paths <- file.path("Rdata", "epred_draws", paste0(epred_ids, ".rds"))
current_epred_md5 <- stats::setNames(unname(tools::md5sum(epred_paths)), epred_ids)
if (!identical(manifest$epred_cache_md5[epred_ids], current_epred_md5)) {
  stop("Orchard figures are older than one or more active prediction caches.", call. = FALSE)
}

required_figures <- file.path(
  "Rdata", "figures", "publication", "orchard",
  paste0(epred_ids, "_orchard_combined.png")
)
missing_figures <- required_figures[!file.exists(required_figures)]
if (length(missing_figures)) {
  stop("Missing v2 orchard figures: ", paste(missing_figures, collapse = ", "),
       call. = FALSE)
}

message("Summary artifact check passed: posterior-mean/equal-tail specification v2.")
