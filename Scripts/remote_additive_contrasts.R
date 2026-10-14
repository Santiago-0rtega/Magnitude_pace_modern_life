# Pairwise moderator contrasts and elapsed-time slopes for the time-adjusted
# (additive) location-scale models, <moderator> + <timevar>, fitted by
# Scripts/remote_fit_one_additive.R. Read-only with respect to the fits.
#
# Contrasts are posterior differences of the moderator coefficients, i.e.
# conditional on the elapsed-time covariate (the reference level contributes
# 0). Location contrasts are on the lnM scale; scale contrasts on the
# log(residual SD) scale. Direction: for levels ordered as in the fitted
# factor (reference first), each pair i < j is reported as level_j - level_i.
# Summaries: posterior mean, 95% equal-tailed credible interval (2.5% and
# 97.5% quantiles), probability of direction (pd), over all post-warmup draws.
#
# Output: <output root>/tables/additive_contrasts.csv
# Run from the repository root (output root set by options(pace.out_root)):
#   Rscript Scripts/remote_additive_contrasts.R

suppressPackageStartupMessages({
  library(here)
  library(posterior)
})
source(here::here("Scripts", "01_paths.R"))

moderators <- c("disturbance", "design", "trait_type", "genphen")
timevars   <- c("log10_years", "log10_generations")

summ <- function(x) {
  c(estimate  = mean(x),
    lower.CrI = unname(stats::quantile(x, 0.025)),
    upper.CrI = unname(stats::quantile(x, 0.975)),
    pd        = max(mean(x > 0), mean(x < 0)))
}

coef_name <- function(moderator, level) {
  paste0(moderator, gsub("[^[:alnum:]_.]", "", level))
}

rows <- list()
for (moderator in moderators) for (timevar in timevars) {
  label <- paste0(moderator, "_plus_", timevar)
  path  <- dir_out("models", paste0(label, "_ls_additive.rds"))
  if (!file.exists(path)) stop("Missing fit: ", label)
  fit   <- readRDS(path)
  levs  <- levels(factor(fit$data[[moderator]]))
  draws <- posterior::as_draws_df(fit)
  n_draws <- nrow(draws)

  get_b <- function(component, level) {
    if (identical(level, levs[[1]])) return(rep(0, n_draws))
    v <- paste0(if (component == "scale") "b_sigma_" else "b_",
                coef_name(moderator, level))
    if (!v %in% names(draws)) stop("Coefficient not found: ", v, " in ", label)
    as.numeric(draws[[v]])
  }

  meta <- data.frame(
    model_id  = label, moderator = moderator, timevar = timevar,
    k         = nrow(fit$data),
    n_refs    = length(unique(as.character(fit$data$ref_id))),
    n_species = length(unique(as.character(fit$data$sp_ncbi))),
    n_draws   = n_draws, stringsAsFactors = FALSE
  )

  for (component in c("location", "scale")) {
    v <- paste0(if (component == "scale") "b_sigma_" else "b_", timevar)
    if (!v %in% names(draws)) stop("Slope not found: ", v, " in ", label)
    s <- summ(as.numeric(draws[[v]]))
    rows[[length(rows) + 1]] <- cbind(meta, data.frame(
      component = component, type = "slope", level_a = NA_character_,
      level_b = NA_character_, contrast = timevar, t(s),
      stringsAsFactors = FALSE))

    for (i in seq_len(length(levs) - 1L)) for (j in (i + 1L):length(levs)) {
      d <- get_b(component, levs[[j]]) - get_b(component, levs[[i]])
      s <- summ(d)
      rows[[length(rows) + 1]] <- cbind(meta, data.frame(
        component = component, type = "contrast",
        level_a = levs[[j]], level_b = levs[[i]],
        contrast = paste(levs[[j]], "-", levs[[i]]), t(s),
        stringsAsFactors = FALSE))
    }
  }
  cat("Done:", label, "| k =", meta$k, "| draws =", n_draws,
      "| levels:", paste(levs, collapse = ", "), "\n")
  rm(fit, draws); invisible(gc())
}

out <- do.call(rbind, rows)
rownames(out) <- NULL
out_path <- dir_out("tables", "additive_contrasts.csv")
utils::write.csv(out, out_path, row.names = FALSE)
cat("Wrote", nrow(out), "rows to tables/additive_contrasts.csv\n")
