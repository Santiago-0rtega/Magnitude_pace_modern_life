setwd(here::here())
library(here)
source(here::here("chapters", "Rscripts", "00_packages_fit.R"))
source(here::here("chapters", "Rscripts", "01_paths.R"))
source(here::here("chapters", "Rscripts", "05_phylogeny.R"))
source(here::here("chapters", "Rscripts", "07_model_formulas.R"))
source(here::here("chapters", "Rscripts", "08_fit_or_read_model.R"))

refit_models <- FALSE

dat_es   <- readRDS(dir_out("effect_sizes", "proceed_lnm_safe.rds"))
A_full   <- readRDS(dir_out("phylogeny", "proceed_A_matrix.rds"))
name_map <- readRDS(dir_out("phylogeny", "proceed_name_map.rds"))

if (!is.null(name_map)) {
  dat_es <- apply_phylo_name_map(dat_es, name_map)
} else {
  dat_es$sp_ncbi_canonical <- dat_es$sp_ncbi
}

# --- Prepare model data ---
dat_model <- dat_es[!is.na(dat_es$disturbance) & !is.na(dat_es$log10_years), ]
dat_model <- dat_model[!is.na(dat_model$sp_ncbi_canonical), ]
dat_model$es_id_model <- factor(dat_model$es_id_model)

cat("n =", nrow(dat_model), "\n")
cat("disturbance levels:\n"); print(table(dat_model$disturbance))

phylo_prep <- prepare_phylo_and_data(dat_model, A_full, label = "m12")
dat_model  <- phylo_prep$dat_model
A_mod      <- phylo_prep$A

V <- metafor::vcalc(vi_lnM_safe, cluster = ref_id,
                    obs = es_id_model, data = dat_model, rho = 0.5)
rownames(V) <- colnames(V) <- dat_model$es_id_model

# --- Interaction formula ---
formula_m12 <- brms::bf(
  yi_lnM_safe ~ log10_years * disturbance +
    (1 | ref_id) +
    (1 | gr(sp_ncbi_canonical, cov = A)) +
    (1 | gr(es_id_model, cov = V)),
  sigma ~ log10_years * disturbance
)

priors_m12 <- c(
  brms::prior(normal(0, 1),   class = b),
  brms::prior(normal(0, 1),   class = b, dpar = sigma),
  brms::prior(normal(0, 0.5), class = sd),
  brms::prior(constant(1),    class = sd, group = es_id_model)
)

fit_m12 <- fit_or_read_model(
  model_name = "m12_ls_years_x_disturbance",
  fit_fun = function() {
    brms::brm(
      formula = formula_m12,
      data    = dat_model,
      data2   = list(A = A_mod, V = V),
      prior   = priors_m12,
      chains  = 4, iter = 4000, warmup = 2000, cores = 4,
      backend = "cmdstanr",
      control = list(adapt_delta = 0.97, max_treedepth = 15)
    )
  },
  refit = refit_models
)

cat("\n=== m12 summary ===\n")
print(summary(fit_m12))
cat("\nDone.\n")
