setwd(here::here())
source(here::here("Scripts", "01_paths.R"))
library(here)
suppressPackageStartupMessages({ source(here::here("Scripts", "00_packages_fit.R")) })
library(emmeans)

out_dir <- dir_out("summaries")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sink_path <- file.path(out_dir, "model_summaries_emmeans.txt")
sink(sink_path, split = TRUE)

# ── Helper: sigma marginals + contrasts from posterior draws ──────────────────
sigma_emmeans <- function(fit, moderator) {
  fe       <- brms::fixef(fit)
  draws    <- brms::as_draws_df(fit)
  ref_lev  <- levels(factor(fit$data[[moderator]]))[1]

  # All levels
  all_levs <- unique(as.character(fit$data[[moderator]]))
  all_levs <- all_levs[!is.na(all_levs)]

  sig_int_draws <- draws[["b_sigma_Intercept"]]

  lev_draws <- list()
  lev_draws[[ref_lev]] <- sig_int_draws

  for (lv in all_levs[all_levs != ref_lev]) {
    col <- paste0("b_sigma_", moderator, lv)
    if (col %in% names(draws)) {
      lev_draws[[lv]] <- sig_int_draws + draws[[col]]
    }
  }

  # Marginal summary
  summ <- do.call(rbind, lapply(names(lev_draws), function(lv) {
    d <- lev_draws[[lv]]
    data.frame(level = lv, emmean = mean(d),
               lower.HPD = coda::HPDinterval(coda::as.mcmc(d))[1],
               upper.HPD = coda::HPDinterval(coda::as.mcmc(d))[2])
  }))

  # All pairwise contrasts
  levs <- names(lev_draws)
  pairs <- do.call(rbind, lapply(seq_along(levs), function(i) {
    do.call(rbind, lapply(seq_along(levs), function(j) {
      if (j <= i) return(NULL)
      d <- lev_draws[[levs[i]]] - lev_draws[[levs[j]]]
      data.frame(
        contrast   = paste(levs[i], "-", levs[j]),
        estimate   = mean(d),
        lower.HPD  = coda::HPDinterval(coda::as.mcmc(d))[1],
        upper.HPD  = coda::HPDinterval(coda::as.mcmc(d))[2],
        p_nonzero  = mean(d > 0)
      )
    }))
  }))

  list(marginals = summ, contrasts = pairs)
}

# ── Model list ────────────────────────────────────────────────────────────────
models <- list(
  list(id="m00", file="m00_ls_intercept_only",    moderator=NULL,                type="intercept"),
  list(id="m01", file="m01_ls_disturbance",        moderator="disturbance",       type="categorical"),
  list(id="m02", file="m02_ls_design",             moderator="design",            type="categorical"),
  list(id="m03", file="m03_ls_log10_years",        moderator="log10_years",       type="continuous"),
  list(id="m04", file="m04_ls_log10_generations",  moderator="log10_generations", type="continuous"),
  list(id="m05", file="m05_ls_trait_type",         moderator="trait_type",        type="categorical"),
  list(id="m06", file="m06_ls_taxa",               moderator="taxa",              type="categorical"),
  list(id="m07", file="m07_ls_genphen",            moderator="genphen",           type="categorical"),
  list(id="m08", file="m08_ls_env_change",         moderator="env_change",        type="categorical"),
  list(id="m09", file="m09_ls_data_type",          moderator="data_type",         type="categorical"),
  list(id="m10", file="m10_ls_transf_data",        moderator="transf_data",       type="categorical"),
  list(id="m11", file="m11_ls_data_scale",         moderator="data_scale",        type="categorical")
)

# ── Main loop ─────────────────────────────────────────────────────────────────
for (m in models) {
  rds_path <- dir_out("models", paste0(m$file, ".rds"))
  if (!file.exists(rds_path)) {
    cat("\n\n========================================\n")
    cat(m$id, "-- SKIPPED (RDS not found)\n")
    cat("========================================\n")
    next
  }

  cat("\n\n========================================\n")
  cat(m$id, "(", m$file, ")\n")
  cat("========================================\n")

  fit <- readRDS(rds_path)

  # Convergence flags
  rhats   <- brms::rhat(fit)
  max_rhat <- max(rhats, na.rm = TRUE)
  n_chains <- fit$fit@sim$chains
  cat(sprintf("Chains: %d | Max Rhat: %.3f\n", n_chains, max_rhat))
  if (max_rhat > 1.01) cat("  ** CONVERGENCE WARNING: Rhat > 1.01 **\n")

  cat("\n--- brms summary ---\n")
  print(summary(fit))

  cat("\n--- Fixed effects (location + scale) ---\n")
  print(round(brms::fixef(fit), 4))

  cat("\n--- Random effects ---\n")
  print(brms::VarCorr(fit))

  # Overall intercept stats
  draws_int <- brms::as_draws_df(fit)[["b_Intercept"]]
  cat(sprintf("\nLocation intercept: mean=%.4f, 95%% CrI [%.4f, %.4f], P(>0)=%.3f\n",
              mean(draws_int), quantile(draws_int, 0.025), quantile(draws_int, 0.975),
              mean(draws_int > 0)))

  if (m$type == "intercept") {
    draws_sig <- brms::as_draws_df(fit)[["b_sigma_Intercept"]]
    cat(sprintf("Sigma intercept:    mean=%.4f, 95%% CrI [%.4f, %.4f]\n",
                mean(draws_sig), quantile(draws_sig, 0.025), quantile(draws_sig, 0.975)))
    cat(sprintf("exp(sigma):         mean=%.4f (residual SD)\n", mean(exp(draws_sig))))
  }

  if (m$type == "categorical") {
    # Location emmeans
    cat("\n--- Location: marginal means (emmeans) ---\n")
    em_loc <- tryCatch(
      emmeans::emmeans(fit, as.formula(paste("~", m$moderator)), epred = TRUE, re_formula = NA),
      error = function(e) { cat("emmeans failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(em_loc)) {
      print(summary(em_loc))
      cat("\n--- Location: pairwise contrasts ---\n")
      pw <- emmeans::contrast(em_loc, method = "pairwise", adjust = "none")
      print(summary(pw, infer = TRUE))
    }

    # Sigma emmeans — draws-based (correct for location-scale brms)
    cat("\n--- Scale (sigma): marginal means from posterior draws ---\n")
    sig_res <- tryCatch(sigma_emmeans(fit, m$moderator),
                        error = function(e) { cat("sigma draws failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(sig_res)) {
      cat("(all values on log-sigma scale; exp() gives residual SD per group)\n")
      print(sig_res$marginals)
      cat("\n--- Scale (sigma): pairwise contrasts ---\n")
      print(sig_res$contrasts)
    }
  }

  if (m$type == "continuous") {
    for (par in c("location", "sigma")) {
      col <- if (par == "location") paste0("b_", m$moderator) else paste0("b_sigma_", m$moderator)
      d   <- tryCatch(brms::as_draws_df(fit)[[col]], error = function(e) NULL)
      if (!is.null(d)) {
        cat(sprintf("\n%s slope: mean=%.4f, 95%% CrI [%.4f, %.4f], P(>0)=%.3f\n",
                    par, mean(d), quantile(d, 0.025), quantile(d, 0.975), mean(d > 0)))
      }
    }
  }
}

sink()
cat("\nSummaries written to:", sink_path, "\n")
