# Fast input validation for Figure 2; run from the repository root.
path <- file.path("data", "bednet_distribution_by_province_2008_2022.xlsx")
if (!file.exists(path)) stop("Missing Figure 2 input: ", path)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Please run source('install_packages.R') first (requires 'readxl').")
}
d <- readxl::read_excel(path, sheet = 1)
stopifnot(identical(names(d)[1], "province"),
          nrow(d) == 26L,
          !anyNA(d$province),
          !anyDuplicated(d$province))
years <- as.character(2009:2024)
missing <- setdiff(years, names(d))
if (length(missing)) stop("Missing year columns: ", paste(missing, collapse = ", "))
values <- unlist(d[years], use.names = FALSE)
if (anyNA(values) || !all(values %in% c(0, 1))) {
  stop("Figure 2 campaign flags must be complete and binary (0/1).")
}
message("Figure 2 input valid: ", nrow(d), " provinces, ", length(years),
        " years, ", sum(values), " flagged province-year campaigns.")
