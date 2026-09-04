# Fitting toolchain. Sourced by the run_*/remote_fit_* scripts instead of
# 00_packages.R, which is deliberately render-only: a book build reads cached
# summaries and must not require brms/cmdstanr merely to print a table.
#
# This file layers the fitting dependencies on top of the render set, so a
# fitting script gets exactly one source() line for packages.

source(here::here("chapters", "Rscripts", "00_packages.R"))

fit_packages <- c("bayesplot", "brms", "cmdstanr", "coda", "emmeans",
                  "posterior", "tidybayes")

missing_fit_packages <- fit_packages[
  !vapply(fit_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_fit_packages) > 0) {
  stop(
    "Missing model-fitting packages: ",
    paste(missing_fit_packages, collapse = ", "),
    ".\ncmdstanr is not on CRAN; install it with:\n",
    '  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))',
    call. = FALSE
  )
}

library(brms)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(bayesplot)

# brms needs a compiled CmdStan, which is a separate install from the R package.
# Fail here with an actionable message rather than deep inside the first brm().
cmdstan_ok <- !inherits(try(cmdstanr::cmdstan_path(), silent = TRUE), "try-error")

if (!cmdstan_ok) {
  stop(
    "CmdStan is not installed or not found.\n",
    "Install it once per machine with:\n",
    "  cmdstanr::check_cmdstan_toolchain(fix = TRUE)\n",
    "  cmdstanr::install_cmdstan(cores = parallel::detectCores())\n",
    "If your home directory is synced (OneDrive/Dropbox), install outside it:\n",
    '  cmdstanr::install_cmdstan(dir = "C:/cmdstan")\n',
    '  then set options(cmdstanr_write_stan_file_dir = ...) / CMDSTAN env var.',
    call. = FALSE
  )
}

options(brms.backend = "cmdstanr")
