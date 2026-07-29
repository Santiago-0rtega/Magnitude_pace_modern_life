# One list is used for setup checks and for the reproducibility chapter.
# Keep it in alphabetical order. Base and recommended R packages are not listed.
project_packages <- c(
  "ape", "bayesplot", "brms", "cli", "cmdstanr", "coda", "dplyr",
  "emmeans", "ggbeeswarm", "ggplot2", "here", "httr2", "janitor",
  "kableExtra", "knitr", "orchaRd", "patchwork", "png", "posterior", "prepR4pcm",
  "purrr", "quarto", "readr", "remotes", "rotl", "scales", "sessioninfo", "stringr",
  "tibble", "tidybayes", "tidyr", "tidyverse"
)

analysis_packages <- setdiff(project_packages, c("quarto", "remotes"))

missing_packages <- analysis_packages[
  !vapply(analysis_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Follow the package installation step in the Reproducibility chapter.",
    call. = FALSE
  )
}

library(tidyverse)
library(here)
library(janitor)
library(brms)
library(cmdstanr)
library(orchaRd)
library(tidybayes)
library(posterior)
library(bayesplot)
library(ape)
library(rotl)
library(prepR4pcm)
library(patchwork)
library(ggbeeswarm)
library(httr2)
library(sessioninfo)
