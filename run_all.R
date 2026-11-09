# Run from the repository root in RStudio: source("run_all.R")
# Final figures: set environment variable FAST_MODE=FALSE (default).
# FAST_MODE=TRUE runs fewer numerical uncertainty draws for preliminary checks.

expected <- file.path("scripts", "01_simulate_main.R")
if (!file.exists(expected)) {
  stop("Run source('run_all.R') from the repository root directory.")
}
for (dir in c("output", "figures", "cache")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

campaign_file <- file.path("data", "bednet_distribution_by_province_2008_2022.xlsx")
if (file.exists(campaign_file)) {
  source(file.path("scripts", "00_check_figure2_input.R"), encoding = "UTF-8")
}
source(file.path("scripts", "01_simulate_main.R"), encoding = "UTF-8")
source(file.path("scripts", "02_figures_main.R"), encoding = "UTF-8")
source(file.path("scripts", "03_figureS1.R"), encoding = "UTF-8")

if (file.exists(campaign_file)) {
  source(file.path("scripts", "04_figure2_provinces.R"), encoding = "UTF-8")
} else {
  warning("Figure 2 skipped because the original 26-province campaign Excel file ",
          "is absent: ", campaign_file, call. = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path("output", "R_sessionInfo.txt"))
message("Finished. Check figures/ and output/.")
