# precompute_heterogeneity.R
# Run ONCE on remote-server (m00 is local + fast there). Computes the location-scale
# heterogeneity decomposition (I^2) for the intercept-only baseline model (m00)
# following Nakagawa et al.'s location-scale meta-analysis framework
# (https://itchyshin.github.io/location-scale_meta-analysis/), extended with the
# phylogenetic variance component present in our models.
#
# Our m00:
#   yi_lnM_safe ~ 1 + (1 | ref_id) + (1 | gr(sp_ncbi, cov = A))
#                   + (1 | gr(es_id_model, cov = V))
#   sigma ~ 1
# Variance components:
#   sigma2_u     = between-study variance (ref_id)                      [location]
#   sigma2_phylo = phylogenetic variance (sp_ncbi, cov = A)             [location]
#   sigma2bar_e  = residual / within-study heterogeneity variance.
#                  Because the scale part carries no study random effect,
#                  Eq. 25 reduces to exp(2 * sigma_Intercept).
#   Vbar         = "typical" sampling variance (Eq. from the tutorial),
#                  computed from the known sampling-variance matrix V.
#
# I^2_x = 100 * sigma2_x / (sigma2_u + sigma2_phylo + sigma2bar_e + Vbar)
#
# `outputs/` is the compute-server workspace. Copy these completed files to the
# matching `Rdata/tables/` and `Rdata/summaries/` book folders before render.
# Output (small): outputs/tables/heterogeneity_m00.csv
#                 outputs/summaries/heterogeneity_m00.rds
#
#   Rscript precompute_heterogeneity.R

setwd(here::here())
source(here::here("Scripts", "01_paths.R"))
suppressMessages({ library(brms); library(dplyr); library(tibble); library(readr) })

model_dir  <- dir_out("models")
tables_dir <- dir_out("tables")
summ_dir   <- dir_out("summaries")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(summ_dir,   showWarnings = FALSE, recursive = TRUE)

m00_path <- file.path(model_dir, "m00_ls_intercept_only.rds")
stopifnot(file.exists(m00_path))
fit <- readRDS(m00_path)

# ── "Typical" sampling variance (Vbar) ───────────────────────────────────────
# Faithful to the tutorial's compute_Vbar(), but (a) reads the sampling-variance
# matrix from data2$V (our name; the tutorial uses "vcv"), and (b) exploits the
# fact that V is diagonal to avoid inverting / forming any n x n dense matrix.
#   Vbar = (n - p_loc) / trace(P),   P = W - W X (X'WX)^{-1} X' W,   W = V^{-1}
# with trace(P) = sum(diag(W)) - trace((X'WX)^{-1} (X'W^2 X)).
compute_Vbar_diag <- function(model) {
  d2  <- model$data2
  Vnm <- names(d2)[vapply(d2, function(m)
           is.matrix(m) && nrow(m) == nobs(model), logical(1))]
  stopifnot(length(Vnm) == 1L)
  v <- diag(d2[[Vnm]])
  w <- 1 / v                                   # diagonal of W = V^{-1}
  X <- brms::make_standata(formula(model), data = model$data,
                           data2 = model$data2)$X
  XtWX  <- crossprod(X * w, X)                  # X' W  X
  XtW2X <- crossprod(X * w)                     # X' W^2 X
  M <- solve(XtWX)
  trace_P <- sum(w) - sum(diag(M %*% XtW2X))
  p_loc <- nrow(brms::fixef(model)) / 2         # equal # loc & scale fixed eff.
  (nobs(model) - p_loc) / trace_P
}

Vbar <- compute_Vbar_diag(fit)

# ── Per-draw variance components → posterior distribution of I^2 ──────────────
draws <- brms::as_draws_df(fit)
sd_ref   <- draws[["sd_ref_id__Intercept"]]
sd_phylo <- draws[["sd_sp_ncbi__Intercept"]]
sig_int  <- draws[["b_sigma_Intercept"]]
stopifnot(!is.null(sd_ref), !is.null(sd_phylo), !is.null(sig_int))

sigma2_u     <- sd_ref^2
sigma2_phylo <- sd_phylo^2
sigma2bar_e  <- exp(2 * sig_int)               # scale part has no study RE
den          <- sigma2_u + sigma2_phylo + sigma2bar_e + Vbar

I2 <- list(
  I2_total = 100 * (sigma2_u + sigma2_phylo + sigma2bar_e) / den,
  I2_study = 100 *  sigma2_u                                / den,  # between-study (ref_id)
  I2_phylo = 100 *  sigma2_phylo                            / den,  # phylogenetic
  I2_resid = 100 *  sigma2bar_e                             / den   # within-study / residual
)

qsum <- function(x) c(median = median(x),
                      l95 = quantile(x, 0.025, names = FALSE),
                      u95 = quantile(x, 0.975, names = FALSE))

het_tbl <- tibble(
  component = c("Total", "Between-study (reference)", "Phylogeny", "Within-study (residual)"),
  symbol    = c("I2_total", "I2_between", "I2_phylo", "I2_within"),
  bind_rows(lapply(I2, qsum))
) |>
  mutate(across(c(median, l95, u95), ~ round(.x, 1)))

meta <- tibble(
  quantity = c("Vbar (typical sampling variance)",
               "sigma2_u (between-study, ref_id)",
               "sigma2_phylo (phylogeny, sp_ncbi)",
               "sigma2bar_e (within-study residual)"),
  value = round(c(Vbar, median(sigma2_u), median(sigma2_phylo),
                  median(sigma2bar_e)), 4)
)

cat("\n=== m00 location-scale heterogeneity (I^2, %) ===\n")
print(as.data.frame(het_tbl))
cat("\n--- variance components (posterior medians) ---\n")
print(as.data.frame(meta))

write_csv(het_tbl, file.path(tables_dir, "heterogeneity_m00.csv"))
saveRDS(list(I2 = het_tbl, components = meta, Vbar = Vbar,
             built = Sys.time()),
        file.path(summ_dir, "heterogeneity_m00.rds"))

cat("\nSaved:", dir_out("tables", "heterogeneity_m00.csv"), "\n")
