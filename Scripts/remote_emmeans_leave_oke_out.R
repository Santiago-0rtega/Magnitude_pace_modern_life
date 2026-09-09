# remote_emmeans_leave_oke_out.R
# Computes location (emmeans) and scale (sigma posterior draws) marginals +
# pairwise contrasts for the leave-Oke-out m01 refit, using the SAME
# functions/approach as precompute_emmeans_contrasts.R (used for every
# primary categorical model). Run on totoro:
#
#   Rscript remote_emmeans_leave_oke_out.R

setwd("/home/ortegara/Documents/PACE")
suppressMessages({
  library(brms); library(emmeans); library(dplyr); library(tibble); library(readr)
})

moderator <- "disturbance"
fit <- readRDS(here::here("outputs", "models", "sensitivity",
                          "m01_ls_disturbance_leave_oke_out.rds"))

out_dir <- here::here("outputs", "tables", "sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Location marginals + pairwise contrasts (identical to precompute_emmeans_contrasts.R) ──
location_contrasts <- function(fit, moderator) {
  em <- emmeans::emmeans(fit, as.formula(paste("~", moderator)),
                         epred = TRUE, re_formula = NA)
  em_draws <- as.matrix(emmeans::as.mcmc.emmGrid(em, names = FALSE))
  em_levels <- as.character(as.data.frame(em)[[1]])
  emm_df <- bind_rows(lapply(seq_along(em_levels), function(i) {
    d <- em_draws[, i]
    tibble(
      level = em_levels[i], emmean = mean(d),
      lower.CrI = unname(quantile(d, 0.025)),
      upper.CrI = unname(quantile(d, 0.975))
    )
  }))

  pw <- emmeans::contrast(em, method = "pairwise", adjust = "none")
  pw_draws <- as.matrix(emmeans::as.mcmc.emmGrid(pw, names = FALSE))
  pw_labels <- as.character(as.data.frame(pw)$contrast)
  ctr_df <- bind_rows(lapply(seq_along(pw_labels), function(i) {
    d <- pw_draws[, i]; p <- mean(d > 0)
    tibble(
      contrast = pw_labels[i], estimate = mean(d),
      lower.CrI = unname(quantile(d, 0.025)),
      upper.CrI = unname(quantile(d, 0.975)), pd = max(p, 1 - p)
    )
  }))
  list(emmeans = emm_df, contrasts = ctr_df)
}

# ── Scale (sigma) marginals + pairwise contrasts, from posterior draws (identical) ──
sigma_contrasts <- function(fit, moderator) {
  draws   <- brms::as_draws_df(fit)
  ref_lev <- levels(factor(fit$data[[moderator]]))[1]
  all_lev <- unique(as.character(fit$data[[moderator]]))
  all_lev <- all_lev[!is.na(all_lev)]

  sig_int <- draws[["b_sigma_Intercept"]]
  lev_draws <- list(); lev_draws[[ref_lev]] <- sig_int
  sig_cols  <- grep(paste0("^b_sigma_", moderator), names(draws), value = TRUE)
  col_key   <- gsub("[^[:alnum:]_]", "", sub(paste0("^b_sigma_", moderator), "", sig_cols))
  for (lv in setdiff(all_lev, ref_lev)) {
    lv_key <- gsub("[^[:alnum:]_]", "", lv)
    hit    <- sig_cols[col_key == lv_key]
    if (length(hit) == 1L) {
      lev_draws[[lv]] <- sig_int + draws[[hit]]
    } else {
      message("  scale: no sigma coefficient matched level '", lv, "' — dropped")
    }
  }

  cri <- function(d) unname(quantile(d, c(0.025, 0.975)))

  emm <- bind_rows(lapply(names(lev_draws), function(lv) {
    d <- lev_draws[[lv]]; h <- cri(d)
    tibble(level = lv, emmean = mean(d), lower.CrI = h[1], upper.CrI = h[2],
           residual_SD = exp(mean(d)))
  }))

  levs  <- names(lev_draws)
  pairs <- bind_rows(lapply(seq_along(levs), function(i)
    bind_rows(lapply(seq_along(levs), function(j) {
      if (j <= i) return(NULL)
      d <- lev_draws[[levs[i]]] - lev_draws[[levs[j]]]; h <- cri(d)
      pdir <- mean(d > 0)
      tibble(contrast = paste(levs[i], "-", levs[j]),
             estimate = mean(d), lower.CrI = h[1], upper.CrI = h[2],
             pd = max(pdir, 1 - pdir))
    }))))
  list(emmeans = emm, contrasts = pairs)
}

loc <- location_contrasts(fit, moderator)
scl <- sigma_contrasts(fit, moderator)

write_csv(loc$emmeans,   file.path(out_dir, "m01_leave_oke_out_location_emmeans.csv"))
write_csv(loc$contrasts, file.path(out_dir, "m01_leave_oke_out_location_contrasts.csv"))
write_csv(scl$emmeans,   file.path(out_dir, "m01_leave_oke_out_scale_emmeans.csv"))
write_csv(scl$contrasts, file.path(out_dir, "m01_leave_oke_out_scale_contrasts.csv"))

message("DONE: emmeans + contrasts written to ", out_dir)
