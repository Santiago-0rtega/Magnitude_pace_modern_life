fit_or_read_model <- function(model_name,
                              fit_fun,
                              model_dir = here::here("Rdata", "models"),
                              refit     = FALSE) {
  if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)

  model_path <- file.path(model_dir, paste0(model_name, ".rds"))

  if (file.exists(model_path) && !refit) {
    message("Reading cached model: ", model_path)
    return(readRDS(model_path))
  }

  # Read-only book render: never fit on the fly — stop with a clear message.
  if (isTRUE(getOption("pace.read_only", FALSE)) && !refit) {
    stop("[read-only render] Model cache not found: ", model_path,
         "\nThe book will not fit models. Sync the .rds, or set refit = TRUE to fit deliberately.",
         call. = FALSE)
  }

  message("Fitting model: ", model_name)
  fit <- fit_fun()

  saveRDS(fit, model_path)
  message("Saved model: ", model_path)
  fit
}

# Default MCMC settings for main analyses
default_mcmc_args <- list(
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(
    adapt_delta  = 0.97,
    max_treedepth = 15
  )
)

# Quick mode for development (fewer iterations)
quick_mcmc_args <- list(
  chains  = 2,
  cores   = 2,
  iter    = 1000,
  warmup  = 500,
  seed    = 123,
  backend = "cmdstanr",
  control = list(
    adapt_delta  = 0.95,
    max_treedepth = 12
  )
)

# Fit one location-scale model for the moderator grid
fit_ls_model <- function(dat_model, formula, priors,
                         V, A = NULL,
                         mcmc_args = default_mcmc_args) {
  data2 <- list(V = V)
  if (!is.null(A)) data2$A <- A

  do.call(
    brms::brm,
    c(
      list(
        formula = formula,
        data    = dat_model,
        data2   = data2,
        family  = gaussian(),
        prior   = priors
      ),
      mcmc_args
    )
  )
}
