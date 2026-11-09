# Run once to install the packages required by the project.
packages <- c("dplyr", "tidyr", "tibble", "purrr", "ggplot2", "scales",
              "patchwork", "readr", "readxl")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
