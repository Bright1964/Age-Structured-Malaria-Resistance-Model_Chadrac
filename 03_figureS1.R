# Supplementary Figure S1: age-specific cases and deaths averted.
# The unpublished age-prevalence Figure 8 is intentionally omitted.

AGE_ANALYSIS_GROUPS <- list(
  "<5" = c(1L, 2L),
  "5-17" = c(3L, 4L),
  "18-42" = c(5L, 6L),
  "43+" = 7L
)
AGE_ANALYSIS_LABELS <- names(AGE_ANALYSIS_GROUPS)

# Common starting resistance level used only for direct product comparisons in
# Figures 8--9. Set AGE_ANALYSIS_USE_COMMON_qR <- FALSE if you instead want the
# product-specific starting values used in Figures 1--7.
AGE_ANALYSIS_USE_COMMON_qR <- TRUE
AGE_ANALYSIS_COMMON_qR <- PUBLICATION_COMMON_INITIAL_qR
AGE_ANALYSIS_MAIN_COVERAGE <- 0.50

age_analysis_start_qR <- function(product) {
  if (isTRUE(AGE_ANALYSIS_USE_COMMON_qR)) {
    return(AGE_ANALYSIS_COMMON_qR)
  }
  initial_qR_for_product(product)
}

# Aggregate a model-age vector into the four manuscript-facing age groups.
aggregate_age_vector <- function(x) {
  vapply(
    AGE_ANALYSIS_GROUPS,
    function(idx) sum(x[idx], na.rm = TRUE),
    numeric(1)
  )
}

# Run one age-specific arm. This function uses the same step_model(), ITN
# coverage function, parameter draws, and national scaling as the existing code.
simulate_age_arm <- function(product,
                             target,
                             drow,
                             itn = TRUE,
                             qR_start = NULL) {

  cov <- if (itn) {
    build_coverage(target, drow$retention_days, drow$peak_cov_mult)
  } else {
    numeric(SIM_LENGTH)
  }

  scale_pop <- NATIONAL_POP / MODEL_HUMAN_POPULATION
  n_age_groups <- length(AGE_ANALYSIS_GROUPS)

  infected_group <- matrix(
    0,
    nrow = SIM_LENGTH,
    ncol = n_age_groups,
    dimnames = list(NULL, AGE_ANALYSIS_LABELS)
  )
  population_group <- infected_group
  clinical_cases_group <- infected_group
  deaths_group <- infected_group
  qR_v <- numeric(SIM_LENGTH)
  coverage_v <- numeric(SIM_LENGTH)

  if (is.null(qR_start)) qR_start <- age_analysis_start_qR(product)
  st <- reset_state_allele_frequency(SPIN, qR_start)

  for (d in seq_len(SIM_LENGTH)) {
    net <- if (itn) {
      get_itn_parameters(product, d, drow, cov[d])
    } else {
      list(
        coverage = 0,
        d_n0 = 0,
        r_n0 = 0,
        d_g = D_PYR,
        selection_fraction = 0,
        nextgen_share = 0
      )
    }

    st <- step_model(
      state = st,
      coverage = net$coverage,
      d_n0 = net$d_n0,
      r_n0 = net$r_n0,
      d_g = net$d_g,
      selection_fraction = net$selection_fraction,
      day = d,
      include_cost = TRUE,
      tmult = TMULT
    )

    infected_by_model_age <- rowSums(st$Iq)
    population_by_model_age <- st$S + infected_by_model_age + st$R

    infected_group[d, ] <- aggregate_age_vector(infected_by_model_age)
    population_group[d, ] <- aggregate_age_vector(population_by_model_age)
    # Clinical case and death counts use the same national burden calibration
    # factor as the existing burden figures.
    clinical_cases_group[d, ] <-
      BURDEN_CALIBRATION_FACTOR * scale_pop *
      aggregate_age_vector(P_CLIN_AGE * st$new_inf_by_age)

    deaths_group[d, ] <-
      BURDEN_CALIBRATION_FACTOR * scale_pop *
      aggregate_age_vector(st$deaths_by_age)

    qR_v[d] <- st$qR
    coverage_v[d] <- st$coverage
  }

  day <- seq_len(SIM_LENGTH)
  year <- STUDY_START_YEAR + (day - 1L) %/% DAYS_PER_YEAR
  month <- (day - 1L) %/% 30L + 1L
  ydec <- 1 + (day - 1) / DAYS_PER_YEAR

  # Monthly prevalence is calculated from infected person-days / total
  # person-days, which is preferable to an unweighted average of daily rates.
  month_levels <- sort(unique(month))
  monthly <- map_dfr(
    month_levels,
    function(mm) {
      idx <- which(month == mm)
      tibble(
        month = mm,
        ydec = mean(ydec[idx]),
        age_group = AGE_ANALYSIS_LABELS,
        prevalence = 100 * colSums(infected_group[idx, , drop = FALSE]) /
          pmax(colSums(population_group[idx, , drop = FALSE]), 1e-12),
        qR = qR_v[idx[length(idx)]],
        coverage = mean(coverage_v[idx])
      )
    }
  )

  year_levels <- sort(unique(year))
  annual <- map_dfr(
    year_levels,
    function(yy) {
      idx <- which(year == yy)
      tibble(
        year = yy,
        sim_year = yy - STUDY_START_YEAR + 1L,
        age_group = AGE_ANALYSIS_LABELS,
        prevalence = 100 * colSums(infected_group[idx, , drop = FALSE]) /
          pmax(colSums(population_group[idx, , drop = FALSE]), 1e-12),
        cases = colSums(clinical_cases_group[idx, , drop = FALSE]),
        deaths = colSums(deaths_group[idx, , drop = FALSE]),
        qR = qR_v[idx[length(idx)]],
        coverage = mean(coverage_v[idx])
      )
    }
  )

  add_meta_age <- function(df) {
    product_value <- as.character(product)
    target_value <- as.numeric(target)
    draw_value <- as.integer(drow$draw)
    start_qR_value <- as.numeric(qR_start)

    df %>%
      mutate(
        product = product_value,
        coverage_target = target_value,
        coverage_label = paste0(as.integer(100 * target_value), "%"),
        draw = draw_value,
        start_qR = start_qR_value,
        age_group = factor(age_group, levels = AGE_ANALYSIS_LABELS)
      )
  }

  list(
    monthly = add_meta_age(monthly),
    annual = add_meta_age(annual)
  )
}

# -----------------------------------------------------------------------------
# Age-specific intervention simulations using a common resistance baseline
# -----------------------------------------------------------------------------
# Cache supplementary computations across repeat plotting runs, invalidating
# them if the code, parameter draws, or relevant settings change.
AGE_CACHE_FILE <- file.path("cache", "figureS1_age_results.rds")
age_cache_signature <- list(
  model_md5 = unname(tools::md5sum(file.path("scripts", "01_simulate_main.R"))),
  supplement_md5 = unname(tools::md5sum(file.path("scripts", "03_figureS1.R"))),
  draws = DRAWS,
  years = c(STUDY_START_YEAR, STUDY_END_YEAR),
  q_initial = AGE_ANALYSIS_COMMON_qR,
  population = c(MODEL_HUMAN_POPULATION, NATIONAL_POP)
)
age_cache <- if (file.exists(AGE_CACHE_FILE)) {
  tryCatch(readRDS(AGE_CACHE_FILE), error = function(e) NULL)
} else NULL
if (is.list(age_cache) && identical(age_cache$signature, age_cache_signature)) {
  message("Using validated age-specific cache for Figure S1.")
  age_monthly <- age_cache$age_monthly
  age_annual <- age_cache$age_annual
  age_impact_pct <- age_cache$age_impact_pct
} else {
  message("Running additional age-structured intervention simulations ...")

  age_grid <- tidyr::crossing(
    product = PRODUCTS,
    coverage_target = COVERAGE_TARGETS,
    draw = seq_len(N_DRAWS)
  )

  age_results <- pmap(
    .l = list(
      age_grid$product,
      age_grid$coverage_target,
      age_grid$draw
    ),
    .f = function(product, coverage_target, draw_index) {
      simulate_age_arm(
        product = product,
        target = coverage_target,
        drow = DRAWS[draw_index, ],
        itn = TRUE,
        qR_start = age_analysis_start_qR(product)
      )
    }
  )

  age_monthly <- map_dfr(age_results, "monthly")
  age_annual <- map_dfr(age_results, "annual")

  # One common no-ITN comparator per draw, initialized at the same q_R used for
  # the direct product comparisons. Product has no effect when ITNs are absent.
  age_no_itn <- map_dfr(
    seq_len(N_DRAWS),
    function(dr) {
      simulate_age_arm(
        product = "Pyrethroid-only nets",
        target = 0,
        drow = DRAWS[dr, ],
        itn = FALSE,
        qR_start = if (isTRUE(AGE_ANALYSIS_USE_COMMON_qR)) {
          AGE_ANALYSIS_COMMON_qR
        } else {
          INITIAL_qR_PYRETHROID
        }
      )$annual %>%
        select(draw, year, sim_year, age_group, no_cases = cases, no_deaths = deaths)
    }
  )

  # Age-specific cumulative cases and deaths averted relative to no ITNs.
  age_impact <- age_annual %>%
    left_join(
      age_no_itn,
      by = c("draw", "year", "sim_year", "age_group")
    ) %>%
    arrange(product, coverage_target, draw, age_group, year) %>%
    group_by(product, coverage_target, coverage_label, draw, age_group) %>%
    mutate(
      cases_averted = pmax(0, no_cases - cases),
      deaths_averted = pmax(0, no_deaths - deaths),
      cum_cases_averted = cumsum(cases_averted),
      cum_deaths_averted = cumsum(deaths_averted)
    ) %>%
    ungroup()

  age_impact_pct <- age_impact %>%
    group_by(product, coverage_target, coverage_label, draw, age_group) %>%
    summarise(
      cases_pct = 100 * sum(cases_averted) / pmax(sum(no_cases), 1e-9),
      deaths_pct = 100 * sum(deaths_averted) / pmax(sum(no_deaths), 1e-9),
      cases_averted_total = sum(cases_averted),
      deaths_averted_total = sum(deaths_averted),
      .groups = "drop"
    )
  saveRDS(list(signature = age_cache_signature,
               age_monthly = age_monthly,
               age_annual = age_annual,
               age_impact_pct = age_impact_pct), AGE_CACHE_FILE)
}

# Export age-specific source data for Supplementary Figure S1.
write_csv(
  age_monthly,
  file.path(OUT_DIR, "age_specific_monthly_trajectories.csv")
)
write_csv(
  age_annual,
  file.path(OUT_DIR, "age_specific_annual_trajectories.csv")
)
write_csv(
  age_impact_pct,
  file.path(OUT_DIR, "age_specific_cases_deaths_averted.csv")
)


# Four panels A/B cases, C/D deaths, matching manuscript supplement.
figureS1_summ <- age_impact_pct %>%
  group_by(product, coverage_label, coverage_target, age_group) %>%
  summarise(
    cases_m = median(cases_pct, na.rm = TRUE),
    cases_lo = q_lo(cases_pct),
    cases_hi = q_hi(cases_pct),
    deaths_m = median(deaths_pct, na.rm = TRUE),
    deaths_lo = q_lo(deaths_pct),
    deaths_hi = q_hi(deaths_pct),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_label = factor(coverage_label, levels = COVERAGE_LABELS),
    age_group = factor(age_group, levels = AGE_ANALYSIS_LABELS)
  )

make_s1_panel <- function(product_name, outcome = c("cases", "deaths"), show_legend = FALSE) {
  outcome <- match.arg(outcome)
  dat <- figureS1_summ %>% filter(product == product_name)

  y_m <- if (outcome == "cases") "cases_m" else "deaths_m"
  y_lo <- if (outcome == "cases") "cases_lo" else "deaths_lo"
  y_hi <- if (outcome == "cases") "cases_hi" else "deaths_hi"
  y_lab <- if (outcome == "cases") "Cases averted over 12 years (%)" else "Deaths averted over 12 years (%)"

  p <- ggplot(
    dat,
    aes(
      x = age_group,
      y = .data[[y_m]],
      fill = coverage_label,
      group = coverage_label
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.76),
      width = 0.68
    ) +
    geom_errorbar(
      aes(ymin = .data[[y_lo]], ymax = .data[[y_hi]]),
      position = position_dodge(width = 0.76),
      width = 0.18,
      linewidth = PUBLICATION_ERRORBAR_WIDTH
    ) +
    scale_fill_manual(values = cov_col, name = "Effective ITN coverage") +
    labs(
      title = product_name,
      x = "Age group (years)",
      y = y_lab
    ) +
    theme_pub()

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

figureS1_A <- make_s1_panel("Pyrethroid-only nets", "cases", TRUE)
figureS1_B <- make_s1_panel("Next-generation nets", "cases", FALSE)
figureS1_C <- make_s1_panel("Pyrethroid-only nets", "deaths", FALSE)
figureS1_D <- make_s1_panel("Next-generation nets", "deaths", FALSE)

figureS1 <- ((figureS1_A | figureS1_B) / (figureS1_C | figureS1_D)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "",
    subtitle = if (isTRUE(AGE_ANALYSIS_USE_COMMON_qR)) {
      paste0("")
    } else {
      "Products use their original product-specific starting resistance levels"
    },
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_fig(
  figureS1,
  stub = "Figure_S1_age_specific_cases_and_deaths_averted",
  width = PUBLICATION_FIG_WIDTH,
  height = PUBLICATION_FIG_HEIGHT_TWO_ROWS
)
show_plot(figureS1)


readr::write_csv(as.data.frame(figureS1_summ),
  file.path(OUT_DIR, "Figure_S1_plot_data.csv"))
message("Supplementary Figure S1 generated from age-specific simulation outputs.")
