args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("n_se", "n_v"))
  stop("Usage: Rscript remote_fit_one_small_study.R n_se|n_v")
predictor <- args[[1]]

setwd("/home/ortegara/Documents/PACE")
library(here)
library(dplyr)
library(metafor)
source(here::here("R", "05_phylogeny.R"))

dat <- readRDS(here::here("outputs", "effect_sizes", "proceed_lnm_safe.rds"))
A_full <- readRDS(here::here("outputs", "phylogeny", "proceed_A_matrix.rds"))
name_map <- readRDS(here::here("outputs", "phylogeny", "proceed_name_map.rds"))

dat <- apply_phylo_name_map(dat, name_map) |>
  filter(
    is.finite(yi_lnM_safe), is.finite(vi_lnM_safe), vi_lnM_safe > 0,
    is.finite(n1), n1 > 0, is.finite(n2), n2 > 0,
    !is.na(sp_ncbi_canonical), nzchar(sp_ncbi_canonical)
  ) |>
  mutate(
    n0 = (n1 * n2) / (n1 + n2),
    n_se = 1 / sqrt(n0),
    n_v = 1 / n0,
    ref_id = droplevels(factor(ref_id)),
    sp_ncbi = factor(sp_ncbi_canonical),
    es_id_model = factor(seq_len(n()))
  )

A <- expand_A_with_unmatched(A_full, levels(dat$sp_ncbi))
fit <- metafor::rma.mv(
  yi = yi_lnM_safe, V = vi_lnM_safe,
  mods = stats::as.formula(paste("~", predictor)),
  random = list(~ 1 | ref_id, ~ 1 | sp_ncbi, ~ 1 | es_id_model),
  R = list(sp_ncbi = A), Rscale = "cor",
  method = "REML", test = "t", sparse = TRUE, data = dat
)

out_dir <- here::here("outputs", "small_study")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(n = nrow(dat), fit = fit), file.path(out_dir, paste0(predictor, ".rds")))

tab <- data.frame(
  model = if (predictor == "n_se")
    "Small-study slope: 1 / sqrt(n0)" else "Adjustment model: 1 / n0",
  term = rownames(fit$beta), estimate = as.numeric(fit$beta),
  std_error = as.numeric(fit$se), ci_lower = as.numeric(fit$ci.lb),
  ci_upper = as.numeric(fit$ci.ub), p_value = as.numeric(fit$pval),
  row.names = NULL
)
write.csv(tab, file.path(out_dir, paste0(predictor, ".csv")), row.names = FALSE)
cat("CONVERGED:", predictor, "n =", nrow(dat), "\n")
print(tab)
