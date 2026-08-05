# One list is used for setup checks and for the reproducibility chapter.
# Keep it in alphabetical order. Base and recommended R packages are not listed.
project_packages <- c(
  "ape", "bayesplot", "brms", "cli", "cmdstanr", "coda", "dplyr",
  "emmeans", "ggbeeswarm", "ggplot2", "here", "httr2", "janitor",
  "kableExtra", "knitr", "metafor", "orchaRd", "patchwork", "png", "posterior", "prepR4pcm",
  "purrr", "quarto", "readr", "remotes", "rotl", "scales", "sessioninfo", "stringr",
  "tibble", "tidybayes", "tidyr", "tidyverse"
)

# Normal renders consume cached model summaries and do not fit models. Keep the
# render-time dependency check separate from the optional fitting toolchain so
# a book build does not require brms/cmdstanr/tidybayes merely to read caches.
render_packages <- c(
  "ape", "dplyr", "ggbeeswarm", "ggplot2", "here", "httr2", "janitor",
  "kableExtra", "knitr", "metafor", "orchaRd", "patchwork", "png", "prepR4pcm",
  "purrr", "readr", "rotl", "scales", "sessioninfo", "stringr", "tibble",
  "tidyr", "tidyverse"
)

missing_packages <- render_packages[
  !vapply(render_packages, requireNamespace, logical(1), quietly = TRUE)
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
library(orchaRd)
library(ape)
library(rotl)
library(prepR4pcm)
library(patchwork)
library(ggbeeswarm)
library(httr2)
library(sessioninfo)
