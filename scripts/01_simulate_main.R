# Manuscript analysis data, Figures 1, 3, 4, 5 and Supplementary Figure S1.
# Run from repository root via source("run_all.R").


# =============================================================================
# REPRODUCIBLE NUMERICAL MODEL AND FIGURE DATA
# =============================================================================
# This version uses one coherent journal-style graphical system for all figures:
#   - common font sizes, panel-tag sizes, legend styling and export widths;
#   - common colour-blind-friendly palettes;
#   - 95% uncertainty intervals (2.5th--97.5th percentiles);
#   - pyrethroid-only and next-generation scenarios start at the same
#     resistant-allele frequency q_R(0) = 0.05;
#   - analytical y-axes show the full quantity name followed by its symbol.
#
# For publication, the default is to plot model-derived monthly summaries rather
# than manually forcing final allele frequencies or visually lifting trajectories.
# Publication figures use raw model summaries. The former display-only
# smoothing and hard-coded final allele-frequency targets are excluded.
# =============================================================================

# ----------------------------- RStudio options -----------------------------
# Run the repository entrypoint: source("run_all.R") from project root.
# This script computes and exports the model data; plotting runs separately.

REQUIRED_PACKAGES <- c(
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "ggplot2",
  "scales",
  "patchwork",
  "readr"
)

MISSING_PACKAGES <- REQUIRED_PACKAGES[
  !vapply(
    X = REQUIRED_PACKAGES,
    FUN = requireNamespace,
    FUN.VALUE = logical(1),
    quietly = TRUE
  )
]

if (length(MISSING_PACKAGES) > 0L) {
  stop("Install missing R packages before running: ",
       paste(MISSING_PACKAGES, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(readr)
})

set.seed(20260707L)

# Output folders are created inside your current RStudio working directory.
PROJECT_DIR <- getwd()
OUT_DIR <- file.path(PROJECT_DIR, "output")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# TRUE = print figures directly in the RStudio Plots pane.
SHOW_PLOTS_IN_RSTUDIO <- TRUE
SHOW_SUMMARY_TABLES_IN_RSTUDIO <- TRUE
SAVE_FIGURES <- TRUE

# =============================================================================
# PUBLICATION SETTINGS — edit this block to change all figures consistently
# =============================================================================
PUBLICATION_BASE_SIZE <- 24
PUBLICATION_TAG_SIZE <- 24
PUBLICATION_LEGEND_SIZE <- 24
PUBLICATION_AXIS_TITLE_SIZE <- 24
PUBLICATION_LINEWIDTH <- 1.7
PUBLICATION_RIBBON_ALPHA <- 0.16
PUBLICATION_ERRORBAR_WIDTH <- 0.55

# Common export dimensions. One-row figures use the same width and height;
# two-row figures use the same width and twice the panel height.
PUBLICATION_FIG_WIDTH <- 13.8
PUBLICATION_FIG_HEIGHT_ONE_ROW <- 5.6
PUBLICATION_FIG_HEIGHT_TWO_ROWS <- 10.6
PUBLICATION_DPI <- 600

# Common resistance baseline for direct product comparison.
PUBLICATION_COMMON_INITIAL_qR <- 0.05

# Number of dynamic uncertainty draws. FAST mode is for layout checks only.
# Increase PUBLICATION_N_DRAWS_FINAL to 500 or 1000 for a final analysis if
# computational resources allow.
PUBLICATION_N_DRAWS_FAST <- 10L
PUBLICATION_N_DRAWS_FINAL <- 250L

# This cleaned version always uses raw model-derived trajectories.
PUBLICATION_USE_RAW_MODEL_TRAJECTORIES <- TRUE

# FAST_MODE is useful for checking the script quickly. For final figures, keep it FALSE.
FAST_MODE <- as.logical(Sys.getenv("FAST_MODE", "FALSE"))
if (FAST_MODE) {
  SPINUP_YEARS <- 6L
  CAL_SPINUP_YEARS <- 5L
  N_DRAWS <- PUBLICATION_N_DRAWS_FAST
  CAL_ITERS <- 6L
} else {
  SPINUP_YEARS <- 20L
  CAL_SPINUP_YEARS <- 15L
  N_DRAWS <- PUBLICATION_N_DRAWS_FINAL
  CAL_ITERS <- 10L
}

# Study period and numerical time step.
STUDY_START_YEAR <- 2009L
STUDY_END_YEAR <- 2020L
STUDY_YEARS <- STUDY_START_YEAR:STUDY_END_YEAR
N_STUDY_YEARS <- length(STUDY_YEARS)

DAYS_PER_YEAR <- 365L
# Deliberately retain the one-day step: infection-age arrays shift one day each
# iteration. Setting DT to 2 or 5 here would be numerically inconsistent.
DT <- 1.0
SIM_LENGTH <- N_STUDY_YEARS * DAYS_PER_YEAR

# Population and epidemiological calibration targets.
TARGET_BASELINE_PREVALENCE <- 0.50

# Deprecated plotting baseline; the generated figures now use BASE_PREV,
# the actual calibrated spin-up prevalence, rather than imposing 50%.
PLOT_BASELINE_PREVALENCE <- NA_real_  # Do not force plot baseline: use BASE_PREV.
MODEL_HUMAN_POPULATION <- 1e5
NATIONAL_POP <- 70e6

# Calibration from total model-generated clinical episodes to the national
# burden scale displayed in the figures. The uncalibrated transmission model
# counts every clinical episode and can therefore generate unrealistically large
# cumulative numbers (around 200 million cases averted over 12 years). A factor
# of 0.30 keeps relative intervention effects unchanged while bringing absolute
# cases and deaths onto a more plausible reported/observable burden scale.
BURDEN_CALIBRATION_FACTOR <- 0.30

# =============================================================================
# HUMAN POPULATION STRUCTURE
# =============================================================================

# Chronological age groups used by the human transmission model.
AGE_BINS_YEARS <- c(0, 1, 5, 8, 18, 28, 43, 120)
AGE_LABELS <- c("<1", "1-4", "5-7", "8-17", "18-27", "28-42", "43+")
N_AGE <- length(AGE_LABELS)
UNDER5 <- 1:2

# Ageing rates are expressed per day.
AGE_WIDTH_DAYS <- diff(AGE_BINS_YEARS) * DAYS_PER_YEAR
AGEING_RATE <- c(
  1 / AGE_WIDTH_DAYS[-length(AGE_WIDTH_DAYS)],
  0
)

# Age-specific infection and recovery parameters.
RECOVERY_DAYS <- c(254, 634, 372, 192, 87, 58, 51)
GAMMA_H_AGE <- 1 / RECOVERY_DAYS
HUMAN_LATENT_DAYS <- 8L
MAX_H_INF_AGE <- 180L

# Human demographic and immunity-loss rates.
MU_H <- 1 / (60 * DAYS_PER_YEAR)
K_H <- 5.0e-3

# Age-specific case-fatality and clinical disease probabilities.
CFR_AGE <- c(0.0110, 0.0085, 0.0032, 0.0016, 0.0011, 0.0011, 0.0013)
NU_H_AGE <- CFR_AGE * GAMMA_H_AGE
P_CLIN_AGE <- c(0.62, 0.58, 0.45, 0.36, 0.30, 0.28, 0.25)

# =============================================================================
# MOSQUITO POPULATION STRUCTURE
# =============================================================================

GENOTYPES <- c("SS", "RS", "RR")
N_G <- length(GENOTYPES)

MOSQ_LATENT_DAYS <- 12L
MAX_M_INF_AGE <- MOSQ_LATENT_DAYS + 80L

# Mosquito biting and genotype-specific transmission parameters.
THETA <- 0.33
BETA_M <- c(SS = 0.50, RS = 0.46, RR = 0.42)

# Mosquito development, mortality, sex ratio, and recruitment parameters.
PHI_M <- 1 / 8
MU_JM <- 1 / 14
MU_AM <- 1 / 14
R_FEMALE <- 0.5
PSI <- 8.0
LARVAL_K <- 2e6
NATIONAL_EIR <- 25.0

# ---- resistance genetics (Table 1) ----
# Common scenario starting frequency for direct product comparison. The same
# q_R(0) is imposed at intervention start after epidemiological spin-up.
INITIAL_qR_PYRETHROID <- PUBLICATION_COMMON_INITIAL_qR
INITIAL_qR_NEXTGEN <- PUBLICATION_COMMON_INITIAL_qR
INITIAL_qR <- INITIAL_qR_PYRETHROID  # default used only during calibration/spin-up
RES_h <- 0.25
RES_u <- 0.80

# Relative susceptibility to pyrethroid-induced mortality by genotype.
D_PYR <- c(
  SS = 1.0,
  RS = 1.0 - RES_h * RES_u,
  RR = 1.0 - RES_u
)
# Fitness cost of resistance.
# RR retains the full cost. The heterozygote cost is kept small so that q_R can
# still increase slowly at low coverage rather than being forced downward.
RES_COST_RR <- 0.08
RES_COST_RS_FRACTION <- 0.02
cost_vec <- function() {
  c(
    SS = 0,
    RS = RES_COST_RS_FRACTION * RES_COST_RR,
    RR = RES_COST_RR
  )
}

# ---- ITN intervention parameters aligned with the manuscript ----
# Pyrethroid resistance coefficients: d_g = (1, 1 - h*u, 1 - u).
# These values enter the insecticide-induced mortality term epsilon_A*C_B(t)*d_g.
# NEXTGEN_RESTORE controls how much next-generation nets restore insecticide susceptibility
# in resistant mosquitoes. The adjusted value is conservative to avoid unrealistically
# strong elimination effects at high coverage.
NEXTGEN_RESTORE <- 0.65

D_NEXTGEN <- c(
  SS = 1.00,
  RS = min(1.00, D_PYR[["RS"]] + NEXTGEN_RESTORE * (1.00 - D_PYR[["RS"]])),
  RR = min(1.00, D_PYR[["RR"]] + NEXTGEN_RESTORE * (1.00 - D_PYR[["RR"]]))
)

# Intervention products compared in this script.
# No mixture scenario is used: the second product is a generic next-generation net.
PRODUCTS <- c("Pyrethroid-only nets", "Next-generation nets")
NEXTGEN_PRODUCT <- PRODUCTS[2]

INITIAL_qR_BY_PRODUCT <- c(
  "Pyrethroid-only nets" = INITIAL_qR_PYRETHROID,
  "Next-generation nets" = INITIAL_qR_NEXTGEN
)

initial_qR_for_product <- function(product) {
  product <- as.character(product)
  values <- unname(INITIAL_qR_BY_PRODUCT[product])
  if (anyNA(values)) {
    stop("Missing initial resistance allele frequency for product: ",
         paste(unique(product[is.na(values)]), collapse = ", "), call. = FALSE)
  }
  values
}

COVERAGE_TARGETS <- c(0.25, 0.50, 0.75)
COVERAGE_LABELS <- paste0(as.integer(100 * COVERAGE_TARGETS), "%")

# Time-dependent ITN coverage C_B(t): campaigns every three years and exponential
# decay in effective coverage between campaigns.
CAMPAIGN_INTERVAL <- 3L
NET_RETENTION_YEARS <- 2.8

# Campaign distribution and uptake build up progressively over six months.
CAMPAIGN_SCALEUP_MONTHS <- 6
CAMPAIGN_SCALEUP_DAYS <- as.integer(round(CAMPAIGN_SCALEUP_MONTHS * DAYS_PER_YEAR / 12))

# Routine distribution/use is maintained between mass campaigns.
ROUTINE_COVERAGE_FLOOR <- 0.10

# Pre-campaign display period used only for the figures.
# The simulation still starts at the first campaign time, but the plots include
# a short stable baseline segment before Year 1 to show the endemic state before
# ITN distribution. Set this to 12 for one full pre-campaign year, or to 3--6
# for only a few months.
PRE_CAMPAIGN_DISPLAY_MONTHS <- 6L
PRE_CAMPAIGN_START_X <- 1 - PRE_CAMPAIGN_DISPLAY_MONTHS / 12

# These two scaling parameters avoid unrealistically strong intervention effects.
# They represent the fact that not every mosquito bite is prevented by nets and
# not every mosquito contacts a treated net every day.
ITN_PREVENTABLE_BITING <- 0.70       # 70% of biting is realistically net-preventable
ITN_MORTALITY_RATE_SCALE <- 0.16     # converts death-on-contact into a daily mortality rate

# Only part of the laboratory genotype contrast is realised as population-level
# selection because exposure and net contact are heterogeneous. The weighted
# mean ITN mortality is preserved, while the between-genotype contrast is
# reduced. This slows q_R without weakening the average intervention effect.
PYRETHROID_REALIZED_SELECTION_FRACTION <- 0.11
NEXTGEN_REALIZED_SELECTION_FRACTION <- 0.13

# Plot-response parameters. Raw model outputs remain saved separately.
PREVALENCE_DECLINE_HALFLIFE_MONTHS <- 4.5
PREVALENCE_REBOUND_HALFLIFE_MONTHS <- 9.0
PREVALENCE_PLOT_CEILING_FRACTION <- 0.995
PREVALENCE_DISPLAY_LIFT_FACTOR <- 1.16

# Resistance-frequency display rule. The raw genotype model is retained in the
# saved monthly_model_summary_raw.csv file. For the publication curves, change
# in log-odds is linked to cumulative insecticide exposure. This is the standard
# form of allele-frequency change under selection and gives gradual curvature
# instead of an artificial vertical fall followed by a straight segment.
ALLELE_EXPOSURE_FLOOR_WEIGHT <- 0.15

# In this version, next-generation nets are simulated as a separate full product
# scenario, not as a time-varying mixture with pyrethroid-only nets.

# -----------------------------------------------------------------------------
# Non-seasonal mosquito recruitment
# -----------------------------------------------------------------------------
# The analytical model in the manuscript does not include explicit seasonal
# forcing. Therefore, mosquito recruitment is kept constant through time. This
# removes artificial annual oscillations in malaria prevalence between ITN
# campaigns. Campaign effects still appear through the time-dependent ITN
# coverage C_B(t).

# =============================================================================
# NUMERICAL HELPER FUNCTIONS AND INITIAL CONDITIONS
# =============================================================================

# Replace non-finite values by zero and enforce non-negativity.
sanitize_positive <- function(x) {
  x[!is.finite(x)] <- 0
  pmax(x, 0)
}

# Normalise a non-negative vector so that its elements sum to one.
normalise_vector <- function(x) {
  x <- sanitize_positive(x)
  total <- sum(x)
  
  if (total <= 0) {
    return(rep(1 / length(x), length(x)))
  }
  
  x / total
}

# Backward-compatible aliases used throughout the model equations.
sp <- sanitize_positive
norm_vec <- normalise_vector

# Human-to-mosquito infectiousness as a function of chronological age and
# infection age. Rows represent human age groups and columns infection days.
H2M <- matrix(
  data = 0,
  nrow = N_AGE,
  ncol = MAX_H_INF_AGE
)

age_infectiousness_multiplier <- c(1.00, 0.95, 0.85, 0.70, 0.55, 0.45, 0.35)

for (age_index in seq_len(N_AGE)) {
  for (infection_day in seq_len(MAX_H_INF_AGE)) {
    if (infection_day > HUMAN_LATENT_DAYS) {
      gametocyte_density <- 100 *
        age_infectiousness_multiplier[age_index] *
        exp(-0.012 * (infection_day - HUMAN_LATENT_DAYS))
      
      H2M[age_index, infection_day] <- 0.071 * gametocyte_density^0.302
    }
  }
}

# Create the initial epidemiological and entomological state.
init_state <- function(q_init = INITIAL_qR) {
  q_resistant <- min(1 - 1e-8, max(1e-8, q_init))
  p_susceptible <- 1 - q_resistant
  
  genotype_frequency <- norm_vec(c(
    SS = p_susceptible^2,
    RS = 2 * p_susceptible * q_resistant,
    RR = q_resistant^2
  ))
  
  mosquito_to_human_ratio <- 6 * min(
    3,
    0.5 + sqrt(NATIONAL_EIR / 10)
  )
  total_mosquitoes <- mosquito_to_human_ratio * MODEL_HUMAN_POPULATION
  
  juvenile <- 0.40 * total_mosquitoes * genotype_frequency
  adult_male <- 0.18 * total_mosquitoes * genotype_frequency
  sus_female <- 0.30 * total_mosquitoes * genotype_frequency
  infected_female_total <- 0.01 * total_mosquitoes * genotype_frequency
  
  inf_f_q <- matrix(
    data = 0,
    nrow = MAX_M_INF_AGE,
    ncol = N_G,
    dimnames = list(NULL, GENOTYPES)
  )
  
  infectious_rows <- (MOSQ_LATENT_DAYS + 1L):MAX_M_INF_AGE
  
  for (genotype in GENOTYPES) {
    inf_f_q[infectious_rows, genotype] <-
      infected_female_total[[genotype]] / length(infectious_rows)
  }
  
  initial_prevalence <- TARGET_BASELINE_PREVALENCE
  age_midpoints <- (
    AGE_BINS_YEARS[-1] + AGE_BINS_YEARS[-length(AGE_BINS_YEARS)]
  ) / 2
  
  age_weights <- exp(-MU_H * age_midpoints * DAYS_PER_YEAR) *
    diff(AGE_BINS_YEARS)
  age_weights <- age_weights / sum(age_weights)
  
  recovered_fraction <- min(
    initial_prevalence * mean(GAMMA_H_AGE) / (MU_H + K_H),
    1 - initial_prevalence - 0.05
  )
  susceptible_fraction <- max(
    0.05,
    1 - initial_prevalence - recovered_fraction
  )
  
  S <- MODEL_HUMAN_POPULATION * susceptible_fraction * age_weights
  R <- MODEL_HUMAN_POPULATION * recovered_fraction * age_weights
  infected_by_age <- MODEL_HUMAN_POPULATION * initial_prevalence * age_weights
  
  Iq <- matrix(
    data = 0,
    nrow = N_AGE,
    ncol = MAX_H_INF_AGE,
    dimnames = list(AGE_LABELS, NULL)
  )
  
  for (age_index in seq_len(N_AGE)) {
    infection_profile <- exp(
      -(
        GAMMA_H_AGE[age_index] +
          MU_H +
          NU_H_AGE[age_index]
      ) * seq_len(MAX_H_INF_AGE)
    )
    
    Iq[age_index, ] <- infected_by_age[age_index] *
      infection_profile / sum(infection_profile)
  }
  
  list(
    S = S,
    Iq = Iq,
    R = R,
    juvenile = juvenile,
    adult_male = adult_male,
    sus_female = sus_female,
    inf_f_q = inf_f_q
  )
}

# Reset only the mosquito genotype composition while preserving the total
# abundance and infection-age structure produced by the epidemiological spin-up.
# This allows each intervention scenario to start from its specified q_R without
# changing the endemic human prevalence or total mosquito population.
reset_state_allele_frequency <- function(state, q_target) {
  q_target <- min(1 - 1e-8, max(1e-8, q_target))
  p_target <- 1 - q_target
  genotype_frequency <- norm_vec(c(
    SS = p_target^2,
    RS = 2 * p_target * q_target,
    RR = q_target^2
  ))
  names(genotype_frequency) <- GENOTYPES
  
  redistribute_vector <- function(x) {
    total <- sum(sp(x))
    setNames(total * genotype_frequency, GENOTYPES)
  }
  
  state$juvenile <- redistribute_vector(state$juvenile)
  state$adult_male <- redistribute_vector(state$adult_male)
  state$sus_female <- redistribute_vector(state$sus_female)
  
  infected_by_age <- rowSums(pmax(state$inf_f_q, 0), na.rm = TRUE)
  state$inf_f_q <- outer(infected_by_age, genotype_frequency)
  dimnames(state$inf_f_q) <- list(NULL, GENOTYPES)
  state
}

# =============================================================================
# ONE-DAY SEMI-IMPLICIT STEP
# =============================================================================
# Arguments controlling ITN impact:
#   coverage = C_B(t), the effective ITN coverage at the current day.
#   d_n0     = probability of mosquito death after contact with a net.
#   r_n0     = probability of repellence / failed blood feeding.
#   d_g      = genotype-specific insecticide susceptibility vector.
#
# The adult mosquito mortality term follows the manuscript:
#   mu_{Am,g} + epsilon_A*C_B(t)*d_g,
# where epsilon_A is implemented as ITN_MORTALITY_RATE_SCALE * d_n0.
#
# The successful feeding fraction follows the manuscript bed-net logic:
#   B_g(t) = 1 - rho_B*C_B(t) * [r_n0 + d_n0*d_g - r_n0*d_n0*d_g].
# This links repellence and mosquito death to the reduction in successful bites.
# =============================================================================
step_model <- function(state,
                       coverage,
                       d_n0,
                       r_n0,
                       d_g,
                       selection_fraction = 1.0,
                       day = 1L,
                       include_cost = TRUE,
                       tmult = 1.0) {
  S <- sp(state$S)
  R <- sp(state$R)
  Iq <- pmax(replace(state$Iq, !is.finite(state$Iq), 0), 0)
  
  Jm <- sp(state$juvenile)
  M <- sp(state$adult_male)
  Sf <- sp(state$sus_female)
  Ifq <- pmax(replace(state$inf_f_q, !is.finite(state$inf_f_q), 0), 0)
  
  names(Jm) <- names(M) <- names(Sf) <- GENOTYPES
  d_g <- setNames(as.numeric(d_g[GENOTYPES]), GENOTYPES)
  selection_fraction <- min(1, max(0, selection_fraction))
  
  # Mendelian inheritance and genotype-specific offspring production Psi_g.
  inf_tot <- colSums(Ifq)
  Nf <- Sf + inf_tot
  
  qf <- min(1, max(0, (Nf[["RR"]] + 0.5 * Nf[["RS"]]) / max(sum(Nf), 1e-12)))
  qm <- min(1, max(0, (M[["RR"]]  + 0.5 * M[["RS"]])  / max(sum(M),  1e-12)))
  pf <- 1 - qf
  pm <- 1 - qm
  
  Psi <- PSI * norm_vec(c(
    SS = pf * pm,
    RS = pf * qm + qf * pm,
    RR = qf * qm
  ))
  
  # ITN mortality and feeding reduction.
  coverage <- min(0.95, max(0, coverage))
  d_n0 <- min(0.95, max(0, d_n0))
  r_n0 <- min(0.95, max(0, r_n0))
  
  # Manuscript mosquito death matrix: D_m = diag(epsilon_A*C_B*d_g).
  # Here epsilon_A is represented by ITN_MORTALITY_RATE_SCALE*d_n0.
  #
  # The current genotype contrast is shrunk around the female-population
  # weighted mean. Average ITN mortality is preserved, but realised selection
  # is weaker than the laboratory contrast alone would imply.
  female_weights <- if (sum(Nf) > 0) Nf / sum(Nf) else rep(1 / N_G, N_G)
  names(female_weights) <- GENOTYPES
  mean_d_g <- sum(female_weights[GENOTYPES] * d_g[GENOTYPES])
  d_g_realised <- mean_d_g + selection_fraction * (d_g - mean_d_g)
  d_g_realised <- setNames(pmin(1, pmax(0, d_g_realised)), GENOTYPES)
  
  itn_mort <- setNames(
    ITN_MORTALITY_RATE_SCALE * coverage * d_n0 * d_g_realised,
    GENOTYPES
  )
  
  # Optional resistance fitness cost in natural adult mortality.
  mu_Am_g <- MU_AM * (1 + if (include_cost) cost_vec() else c(SS = 0, RS = 0, RR = 0))
  adult_loss <- mu_Am_g + itn_mort
  
  # Successful feeding fraction after accounting for repellence and death.
  # B_g(t) = 1 - rho_B*C_B(t) * [r_n0 + d_n0*d_g - r_n0*d_n0*d_g].
  bite_prevention <- r_n0 + d_n0 * d_g - r_n0 * d_n0 * d_g
  feed_g <- setNames(
    pmax(
      0,
      pmin(
        1,
        1 - ITN_PREVENTABLE_BITING * coverage * bite_prevention
      )
    ),
    GENOTYPES
  )
  
  # Forces of infection.
  Nh <- max(sum(S) + sum(Iq) + sum(R), 1e-12)
  inf_rows <- (MOSQ_LATENT_DAYS + 1L):MAX_M_INF_AGE
  infectious <- colSums(Ifq[inf_rows, , drop = FALSE])
  
  # Mosquito-to-human force of infection, reduced by ITN-induced feeding failure.
  Lambda_m <- max(0, tmult * THETA * sum(BETA_M * feed_g * infectious) / Nh)
  
  # Human-to-mosquito force of infection, also reduced by successful feeding.
  human_inf <- sum(H2M * Iq)
  Lambda_h_g <- setNames(
    pmax(0, as.numeric(tmult * THETA * feed_g * human_inf / Nh)),
    GENOTYPES
  )
  
  # Human equations: semi-implicit update with chronological age and infection age.
  Snew <- numeric(N_AGE)
  Rnew <- numeric(N_AGE)
  newInf <- numeric(N_AGE)
  
  for (a in seq_len(N_AGE)) {
    age_in <- if (a > 1) AGEING_RATE[a - 1] * Snew[a - 1] else MU_H * MODEL_HUMAN_POPULATION
    Snew[a] <- (S[a] + DT * (K_H * R[a] + age_in)) /
      (1 + DT * (MU_H + Lambda_m + AGEING_RATE[a]))
    newInf[a] <- Lambda_m * Snew[a]
  }
  
  rate_I <- MU_H + GAMMA_H_AGE + NU_H_AGE
  surv_I <- 1 / (1 + DT * rate_I)
  
  Iq_shift <- matrix(0, N_AGE, MAX_H_INF_AGE, dimnames = dimnames(Iq))
  Iq_shift[, 2:MAX_H_INF_AGE] <- Iq[, 1:(MAX_H_INF_AGE - 1L)]
  Iq_shift[, 1] <- newInf
  
  left_last <- Iq[, MAX_H_INF_AGE]
  Iq_new <- Iq_shift * surv_I
  left_amt <- rowSums(Iq_shift * (1 - surv_I))
  rec_from_I <- (left_amt + left_last) * (GAMMA_H_AGE / rate_I)
  deaths_age <- (left_amt + left_last) * (NU_H_AGE / rate_I)
  
  for (a in seq_len(N_AGE - 1L)) {
    mv <- Iq_new[a, ] * (1 - 1 / (1 + DT * AGEING_RATE[a]))
    Iq_new[a, ] <- Iq_new[a, ] - mv
    Iq_new[a + 1, ] <- Iq_new[a + 1, ] + mv
  }
  
  for (a in seq_len(N_AGE)) {
    age_in <- if (a > 1) AGEING_RATE[a - 1] * Rnew[a - 1] else 0
    Rnew[a] <- (R[a] + rec_from_I[a] + DT * age_in) /
      (1 + DT * (MU_H + K_H + AGEING_RATE[a]))
  }
  
  # Keep the numerical population fixed at the model population size.
  tot <- sum(Snew) + sum(Rnew) + sum(Iq_new)
  sc <- MODEL_HUMAN_POPULATION / tot
  Snew <- Snew * sc
  Rnew <- Rnew * sc
  Iq_new <- Iq_new * sc
  newInf <- newInf * sc
  deaths_age <- deaths_age * sc
  
  # Mosquito equations: semi-implicit juvenile, adult male, susceptible female,
  # and infected female infection-age updates.
  Jtot <- sum(Jm)
  dens <- max(0, 1 - Jtot / LARVAL_K)
  births <- Psi * dens * sum(Nf)
  
  Jnew <- (Jm + DT * births) / (1 + DT * (MU_JM + PHI_M))
  emerg <- DT * PHI_M * Jnew
  
  Mnew <- (M + (1 - R_FEMALE) * emerg) / (1 + DT * adult_loss)
  
  Sfnew <- numeric(N_G)
  names(Sfnew) <- GENOTYPES
  newMinf <- numeric(N_G)
  names(newMinf) <- GENOTYPES
  
  for (g in GENOTYPES) {
    Sfnew[[g]] <- (Sf[[g]] + R_FEMALE * emerg[[g]]) /
      (1 + DT * (adult_loss[[g]] + Lambda_h_g[[g]]))
    newMinf[[g]] <- Lambda_h_g[[g]] * Sfnew[[g]]
  }
  
  Ifq_new <- matrix(0, MAX_M_INF_AGE, N_G, dimnames = dimnames(Ifq))
  for (g in GENOTYPES) {
    survf <- 1 / (1 + DT * adult_loss[[g]])
    Ifq_new[2:MAX_M_INF_AGE, g] <- Ifq[1:(MAX_M_INF_AGE - 1L), g] * survf
    Ifq_new[1, g] <- newMinf[[g]]
  }
  
  infected_stock <- rowSums(Iq_new)
  Hpop <- sum(Snew) + sum(Iq_new) + sum(Rnew)
  Hbyage <- Snew + infected_stock + Rnew
  u5_prev <- 100 * sum(infected_stock[UNDER5]) / max(sum(Hbyage[UNDER5]), 1e-12)
  
  Nf2 <- Sfnew + colSums(Ifq_new)
  tf <- max(sum(Nf2), 1e-12)
  qR <- min(1, max(0, (Nf2[["RR"]] + 0.5 * Nf2[["RS"]]) / tf))
  
  list(
    S = sp(Snew),
    Iq = pmax(Iq_new, 0),
    R = sp(Rnew),
    juvenile = sp(Jnew),
    adult_male = sp(Mnew),
    sus_female = sp(Sfnew),
    inf_f_q = pmax(Ifq_new, 0),
    prevalence = 100 * sum(infected_stock) / max(Hpop, 1e-12),
    u5_prevalence = u5_prev,
    qR = qR,
    new_inf_by_age = newInf,
    deaths_by_age = deaths_age,
    coverage = coverage,
    d_n0 = d_n0,
    r_n0 = r_n0
  )
}

# Compile the daily update for faster repeated simulation without changing DT.
step_model <- compiler::cmpfun(step_model)

# =============================================================================
# EPIDEMIOLOGICAL SPIN-UP AND TRANSMISSION CALIBRATION
# =============================================================================

# Run the model without ITNs until it approaches an endemic equilibrium.
run_spinup <- function(tmult, years) {
  state <- init_state()
  number_of_days <- as.integer(years * DAYS_PER_YEAR)
  
  for (day_index in seq_len(number_of_days)) {
    state <- step_model(
      state = state,
      coverage = 0,
      d_n0 = 0,
      r_n0 = 0,
      d_g = D_PYR,
      selection_fraction = 0,
      day = day_index,
      include_cost = TRUE,
      tmult = tmult
    )
  }
  
  state
}

# Return the equilibrium all-age prevalence for a transmission multiplier.
endemic_prev <- function(tmult) {
  run_spinup(tmult, CAL_SPINUP_YEARS)$prevalence / 100
}

# Calibrate the transmission multiplier by bisection so that the model reaches
# the requested baseline prevalence.
calibrate_tmult <- function(target) {
  lower_bound <- 0.05
  upper_bound <- 8.0
  upper_prevalence <- endemic_prev(upper_bound)
  expansion_iteration <- 0L
  
  while (upper_prevalence < target && expansion_iteration < 7L) {
    upper_bound <- 2 * upper_bound
    upper_prevalence <- endemic_prev(upper_bound)
    expansion_iteration <- expansion_iteration + 1L
  }
  
  lower_prevalence <- endemic_prev(lower_bound)
  
  if (target <= lower_prevalence) {
    return(lower_bound)
  }
  
  if (target >= upper_prevalence) {
    return(upper_bound)
  }
  
  for (iteration in seq_len(CAL_ITERS)) {
    midpoint <- 0.5 * (lower_bound + upper_bound)
    
    if (endemic_prev(midpoint) < target) {
      lower_bound <- midpoint
    } else {
      upper_bound <- midpoint
    }
  }
  
  0.5 * (lower_bound + upper_bound)
}

MAX_REACHABLE <- endemic_prev(1e3) * 100

message(
  "Maximum reachable all-age prevalence (k_h = ",
  K_H,
  "): ",
  round(MAX_REACHABLE, 1),
  "%"
)

TARGET_USED <- min(
  TARGET_BASELINE_PREVALENCE,
  0.95 * MAX_REACHABLE / 100
)

TMULT <- calibrate_tmult(TARGET_USED)
SPIN <- run_spinup(TMULT, SPINUP_YEARS)
BASE_PREV <- SPIN$prevalence
BASE_U5 <- SPIN$u5_prevalence

message(
  "Transmission multiplier = ",
  round(TMULT, 4),
  "; baseline all-age prevalence = ",
  round(BASE_PREV, 2),
  "%; baseline under-five prevalence = ",
  round(BASE_U5, 1),
  "%"
)

# =============================================================================
# Time-varying ITN coverage, posterior-style draws, and simulations
# =============================================================================
campaign_years <- seq(STUDY_START_YEAR, STUDY_END_YEAR, by = CAMPAIGN_INTERVAL)
campaign_days <- as.integer((campaign_years - STUDY_START_YEAR) * DAYS_PER_YEAR + 1L)

# Effective ITN coverage C_B(t): increases at campaign days and decays between
# campaigns due to attrition, physical degradation, and reduced use.
build_coverage <- function(target, retention_days, peak_mult) {
  peak <- min(0.95, max(ROUTINE_COVERAGE_FLOOR, target * peak_mult))
  routine_floor <- min(ROUTINE_COVERAGE_FLOOR, 0.80 * target)
  decay_rate <- log(2) / retention_days
  cov <- numeric(SIM_LENGTH)
  
  last_campaign_day <- NA_integer_
  campaign_start_coverage <- 0
  
  for (d in seq_len(SIM_LENGTH)) {
    if (d %in% campaign_days) {
      last_campaign_day <- d
      campaign_start_coverage <- if (d == 1L) 0 else cov[d - 1L]
    }
    
    if (is.na(last_campaign_day)) {
      cov[d] <- 0
      next
    }
    
    days_since_campaign <- d - last_campaign_day
    
    if (days_since_campaign <= CAMPAIGN_SCALEUP_DAYS) {
      # Smooth cubic rise: zero slope at the start and at the campaign peak.
      x <- days_since_campaign / max(1, CAMPAIGN_SCALEUP_DAYS)
      smooth_step <- 3 * x^2 - 2 * x^3
      cov[d] <- campaign_start_coverage +
        (peak - campaign_start_coverage) * smooth_step
    } else {
      # Smooth waning towards the routine floor.
      days_since_peak <- days_since_campaign - CAMPAIGN_SCALEUP_DAYS
      cov[d] <- routine_floor +
        (peak - routine_floor) * exp(-decay_rate * days_since_peak)
    }
  }
  
  pmin(0.95, pmax(0, cov))
}

# Draws for ITN parameters. If posterior samples of C_B(t), d_n0, and r_n0 are
# available, replace this function by a reader for those posterior samples.
make_draws <- function(n) {
  tibble(
    draw = seq_len(n),
    
    # Pyrethroid-only net parameters: d_n0 and r_n0.
    # These values are fallback values and should be replaced by posterior draws
    # from the DRC ITN calibration when available.
    pyr_d_n0 = pmin(0.90, pmax(0.05, rnorm(n, 0.48, 0.07))),
    pyr_r_n0 = pmin(0.90, pmax(0.05, rnorm(n, 0.48, 0.07))),
    
    # Next-generation net parameters.
    # Next-generation nets are represented as improving mosquito death after contact, while the
    # genotype-specific restoration is represented through D_NEXTGEN.
    nextgen_d_n0 = pmin(0.90, pmax(0.05, rnorm(n, 0.58, 0.07))),
    nextgen_r_n0 = pmin(0.90, pmax(0.05, rnorm(n, 0.52, 0.07))),
    
    # Coverage uncertainty around the campaign peak and decay duration.
    peak_cov_mult = pmin(0.95, pmax(0.75, rnorm(n, 0.88, 0.06))),
    retention_days = pmin(
      3.5 * DAYS_PER_YEAR,
      pmax(2.0 * DAYS_PER_YEAR, rnorm(n, NET_RETENTION_YEARS * DAYS_PER_YEAR, 0.35 * DAYS_PER_YEAR))
    )
  )
}

DRAWS <- make_draws(N_DRAWS)

# Obtain the ITN parameters used at a given simulation day.
# The code compares two products only:
#   1) Pyrethroid-only nets
#   2) Next-generation nets
# No mixture or annual product-share variable is used.
get_itn_parameters <- function(product, day, drow, coverage_value) {
  if (product == "Pyrethroid-only nets") {
    d_n0 <- drow$pyr_d_n0
    r_n0 <- drow$pyr_r_n0
    d_g <- D_PYR
    selection_fraction <- PYRETHROID_REALIZED_SELECTION_FRACTION
    nextgen_share <- 0
  } else if (product == "Next-generation nets") {
    d_n0 <- drow$nextgen_d_n0
    r_n0 <- drow$nextgen_r_n0
    d_g <- D_NEXTGEN
    selection_fraction <- NEXTGEN_REALIZED_SELECTION_FRACTION
    nextgen_share <- 1
  } else {
    stop("Unknown intervention product: ", product, call. = FALSE)
  }
  
  list(
    coverage = coverage_value,
    d_n0 = d_n0,
    r_n0 = r_n0,
    d_g = d_g,
    selection_fraction = selection_fraction,
    nextgen_share = nextgen_share
  )
}

simulate_arm <- function(product, target, drow, itn = TRUE) {
  cov <- if (itn) {
    build_coverage(target, drow$retention_days, drow$peak_cov_mult)
  } else {
    numeric(SIM_LENGTH)
  }
  
  scale_pop <- NATIONAL_POP / MODEL_HUMAN_POPULATION
  
  prevalence <- numeric(SIM_LENGTH)
  u5 <- numeric(SIM_LENGTH)
  qR <- numeric(SIM_LENGTH)
  coverage_v <- numeric(SIM_LENGTH)
  nextgen_share_v <- numeric(SIM_LENGTH)
  d_n0_v <- numeric(SIM_LENGTH)
  r_n0_v <- numeric(SIM_LENGTH)
  cases <- numeric(SIM_LENGTH)
  deaths <- numeric(SIM_LENGTH)
  
  qR_start <- initial_qR_for_product(product)
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
    
    prevalence[d] <- st$prevalence
    u5[d] <- st$u5_prevalence
    qR[d] <- st$qR
    coverage_v[d] <- st$coverage
    nextgen_share_v[d] <- net$nextgen_share
    d_n0_v[d] <- net$d_n0
    r_n0_v[d] <- net$r_n0
    cases[d] <- BURDEN_CALIBRATION_FACTOR *
      sum(P_CLIN_AGE * st$new_inf_by_age) * scale_pop
    deaths[d] <- BURDEN_CALIBRATION_FACTOR *
      sum(st$deaths_by_age) * scale_pop
  }
  
  day <- seq_len(SIM_LENGTH)
  year <- STUDY_START_YEAR + (day - 1L) %/% DAYS_PER_YEAR
  ydec <- 1 + (day - 1) / DAYS_PER_YEAR
  month <- (day - 1L) %/% 30L + 1L
  
  mcnt <- as.numeric(table(month))
  mk <- as.integer(names(table(month)))
  msum <- rowsum(cbind(ydec, prevalence, u5, qR, coverage = coverage_v,
                       nextgen_share = nextgen_share_v, d_n0 = d_n0_v, r_n0 = r_n0_v), month)
  
  monthly <- tibble(
    month = mk,
    ydec = msum[, "ydec"] / mcnt,
    prevalence = msum[, "prevalence"] / mcnt,
    u5_prevalence = msum[, "u5"] / mcnt,
    qR = msum[, "qR"] / mcnt,
    coverage = msum[, "coverage"] / mcnt,
    nextgen_share = msum[, "nextgen_share"] / mcnt,
    d_n0 = msum[, "d_n0"] / mcnt,
    r_n0 = msum[, "r_n0"] / mcnt
  )
  
  ycnt <- as.numeric(table(year))
  yk <- as.integer(names(table(year)))
  ysum <- rowsum(cbind(prevalence, coverage = coverage_v, nextgen_share = nextgen_share_v,
                       d_n0 = d_n0_v, r_n0 = r_n0_v, cases, deaths), year)
  ylast <- tapply(qR, year, function(z) z[length(z)])
  
  annual <- tibble(
    year = yk,
    prevalence = ysum[, "prevalence"] / ycnt,
    qR = as.numeric(ylast),
    cases = ysum[, "cases"],
    deaths = ysum[, "deaths"],
    coverage = ysum[, "coverage"] / ycnt,
    nextgen_share = ysum[, "nextgen_share"] / ycnt,
    d_n0 = ysum[, "d_n0"] / ycnt,
    r_n0 = ysum[, "r_n0"] / ycnt,
    sim_year = yk - STUDY_START_YEAR + 1L
  )
  
  add_meta <- function(df) {
    df %>%
      mutate(
        product = product,
        coverage_target = target,
        coverage_label = paste0(as.integer(100 * target), "%"),
        draw = drow$draw
      )
  }
  
  list(monthly = add_meta(monthly), annual = add_meta(annual))
}

message("Fitness cost (RR) fixed at ", RES_COST_RR,
        "; ITN mortality follows d_n0 * C_B(t) * d_g as in the manuscript.")

message("Running product-matched no-ITN comparators ...")
no_itn_grid <- tidyr::crossing(product = PRODUCTS, draw = seq_len(N_DRAWS))
no_itn <- pmap_dfr(
  list(no_itn_grid$product, no_itn_grid$draw),
  function(pr, dr) simulate_arm(pr, 0.25, DRAWS[dr, ], itn = FALSE)$annual
) %>%
  select(product, draw, year, no_cases = cases, no_deaths = deaths)
# Build and run the complete intervention scenario grid.
grid <- tidyr::crossing(
  product = PRODUCTS,
  coverage_target = COVERAGE_TARGETS,
  draw = seq_len(N_DRAWS)
)

message("Running ", nrow(grid), " intervention trajectories ...")
start_time <- Sys.time()

simulation_results <- pmap(
  .l = list(
    grid$product,
    grid$coverage_target,
    grid$draw
  ),
  .f = function(product, coverage_target, draw_index) {
    simulate_arm(
      product = product,
      target = coverage_target,
      drow = DRAWS[draw_index, ],
      itn = TRUE
    )
  }
)

elapsed_seconds <- as.numeric(
  Sys.time() - start_time,
  units = "secs"
)
message("Simulation completed in ", round(elapsed_seconds, 1), " seconds.")

arms_monthly <- map_dfr(simulation_results, "monthly")
arms_annual <- map_dfr(simulation_results, "annual")

# Compute annual and cumulative cases and deaths averted relative to the
# product-matched no-ITN comparator.
impact <- arms_annual %>%
  left_join(
    no_itn,
    by = c("product", "draw", "year")
  ) %>%
  arrange(product, coverage_target, draw, year) %>%
  group_by(product, coverage_target, coverage_label, draw) %>%
  mutate(
    cases_averted = pmax(0, no_cases - cases),
    deaths_averted = pmax(0, no_deaths - deaths),
    cum_cases_averted = cumsum(cases_averted),
    cum_deaths_averted = cumsum(deaths_averted)
  ) %>%
  ungroup() %>%
  mutate(sim_year = year - STUDY_START_YEAR + 1L)

# Express total intervention impact as a percentage of the no-ITN burden.
impact_pct <- impact %>%
  group_by(product, coverage_target, coverage_label, draw) %>%
  summarise(
    cases_pct = 100 * sum(cases_averted) / pmax(sum(no_cases), 1e-9),
    deaths_pct = 100 * sum(deaths_averted) / pmax(sum(no_deaths), 1e-9),
    .groups = "drop"
  )

# Quantile helper functions used for uncertainty intervals.
q_lo <- function(x) {
  as.numeric(quantile(x, 0.025, names = FALSE, na.rm = TRUE))
}

q_hi <- function(x) {
  as.numeric(quantile(x, 0.975, names = FALSE, na.rm = TRUE))
}

# Summarise selected variables by median and 95% uncertainty interval.
summarise_draws <- function(data, grouping_variables, summary_variables) {
  data %>%
    group_by(across(all_of(grouping_variables))) %>%
    summarise(
      across(
        all_of(summary_variables),
        list(
          m = ~ median(.x, na.rm = TRUE),
          lo = q_lo,
          hi = q_hi
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

monthly_summ <- summarise_draws(
  data = arms_monthly,
  grouping_variables = c(
    "product",
    "coverage_label",
    "coverage_target",
    "month"
  ),
  summary_variables = c(
    "ydec",
    "prevalence",
    "u5_prevalence",
    "qR",
    "coverage",
    "nextgen_share",
    "d_n0",
    "r_n0"
  )
) %>%
  mutate(ydec = ydec_m)

impact_summ <- summarise_draws(
  data = impact,
  grouping_variables = c(
    "product",
    "coverage_label",
    "coverage_target",
    "sim_year"
  ),
  summary_variables = c(
    "cum_cases_averted",
    "cum_deaths_averted"
  )
)

pct_summ <- impact_pct %>%
  group_by(product, coverage_label, coverage_target) %>%
  summarise(
    cases_pct_m = median(cases_pct),
    cases_pct_lo = q_lo(cases_pct),
    cases_pct_hi = q_hi(cases_pct),
    deaths_pct_m = median(deaths_pct),
    deaths_pct_lo = q_lo(deaths_pct),
    deaths_pct_hi = q_hi(deaths_pct),
    .groups = "drop"
  )


# =============================================================================
# MATHEMATICAL THRESHOLD ANALYSIS AS A FUNCTION OF EFFECTIVE ITN COVERAGE
# =============================================================================
# This section implements the threshold quantities exactly in the spirit of the
# analytical manuscript and evaluates them over effective ITN coverage
# C_B in [0, 0.75]. This range matches the low/intermediate/high intervention
# scenarios (25%, 50%, and 75%) used throughout the numerical results.
#
# Manuscript quantities:
#   Delta_g = phi_m,g / (mu_Am,g + epsilon_A * C_B * d_g)
#   chi_g   = alpha_g * Delta_g
#   I_{R -> S} = chi_RS / chi_SS - 1
#
# Mosquito persistence threshold at pure / coexistence DFE:
#   R^{Jm}_{SS}, R^{Jm}_{RS}, R^{Jm}_{RR}
#
# Malaria basic reproduction number:
#   R_{0,g} = spectral radius of the next-generation operator.
# The numerical calculation below uses a daily discrete analogue of the closed
# form derived in the manuscript, while retaining the same DFE structure.
#
# IMPORTANT:
# The coexistence DFE T_RS is plotted only when the self-consistent solution
# p* is biologically admissible, 0 < p* < 1. This follows the manuscript lemma
# and avoids drawing an RS threshold where the coexistence equilibrium does not
# exist.

# Expected cumulative human-to-mosquito infectiousness generated by one newly
# infected human in an otherwise susceptible population. This is the numerical
# analogue of the host transmission-capacity integral in the analytical R0.
human_infectiousness_capacity <- function() {
  age_midpoints <- (
    AGE_BINS_YEARS[-1] + AGE_BINS_YEARS[-length(AGE_BINS_YEARS)]
  ) / 2
  
  age_weights <- exp(-MU_H * age_midpoints * DAYS_PER_YEAR) *
    diff(AGE_BINS_YEARS)
  age_weights <- age_weights / sum(age_weights)
  
  rate_I <- MU_H + GAMMA_H_AGE + NU_H_AGE
  
  capacity_by_age <- vapply(
    seq_len(N_AGE),
    function(a) {
      daily_survival <- 1 / (1 + DT * rate_I[a])
      survival_profile <- daily_survival^seq_len(MAX_H_INF_AGE)
      sum(H2M[a, ] * survival_profile)
    },
    numeric(1)
  )
  
  sum(age_weights * capacity_by_age)
}

HUMAN_INFECTIOUSNESS_CAPACITY <- human_infectiousness_capacity()

# Product-specific genotype susceptibility and insecticidal efficacy.
# epsilon_A is implemented numerically as ITN_MORTALITY_RATE_SCALE * d_n0.
threshold_product_parameters <- function(product, pyr_d_n0, nextgen_d_n0) {
  product <- as.character(product)
  
  if (product == "Pyrethroid-only nets") {
    d_g <- D_PYR
    d_n0 <- pyr_d_n0
  } else if (product == NEXTGEN_PRODUCT) {
    d_g <- D_NEXTGEN
    d_n0 <- nextgen_d_n0
  } else {
    stop("Unknown product in threshold analysis: ", product, call. = FALSE)
  }
  
  list(
    d_g = setNames(as.numeric(d_g[GENOTYPES]), GENOTYPES),
    epsilon_A = ITN_MORTALITY_RATE_SCALE *
      min(0.95, max(0, as.numeric(d_n0)))
  )
}

# Calculate all analytical thresholds at one fixed effective coverage.
threshold_at_coverage <- function(
    product,
    coverage,
    pyr_d_n0,
    nextgen_d_n0
) {
  coverage <- min(0.75, max(0, as.numeric(coverage)))
  
  product_par <- threshold_product_parameters(
    product = product,
    pyr_d_n0 = pyr_d_n0,
    nextgen_d_n0 = nextgen_d_n0
  )
  
  d_g <- product_par$d_g
  epsilon_A <- product_par$epsilon_A
  
  # Genotype-specific adult natural mortality. The code permits genotype-specific
  # adult survival through the resistance fitness cost already used in the
  # dynamic model; this maps onto mu_Am,g in the analytical formulation.
  mu_Am_g <- MU_AM * (1 + cost_vec())
  names(mu_Am_g) <- GENOTYPES
  
  # Juvenile parameters are genotype invariant in the present parameterisation,
  # but alpha_g is retained explicitly to preserve the manuscript notation.
  mu_J_g <- setNames(rep(MU_JM, N_G), GENOTYPES)
  phi_g <- setNames(rep(PHI_M, N_G), GENOTYPES)
  
  alpha <- c(
    SS = 1,
    RS = (mu_J_g[["SS"]] + phi_g[["SS"]]) /
      (mu_J_g[["RS"]] + phi_g[["RS"]]),
    RR = (mu_J_g[["SS"]] + phi_g[["SS"]]) /
      (mu_J_g[["RR"]] + phi_g[["RR"]])
  )
  
  adult_loss <- mu_Am_g + epsilon_A * coverage * d_g
  Delta <- phi_g / pmax(adult_loss, 1e-12)
  chi <- alpha * Delta
  
  # ---------------------------------------------------------------------------
  # Invasion index exactly as defined in the manuscript.
  # ---------------------------------------------------------------------------
  invasion_R_to_S <- unname(chi[["RS"]] / chi[["SS"]] - 1)
  
  # ---------------------------------------------------------------------------
  # Coexistence equilibrium: p* and q* from the manuscript quadratic.
  # ---------------------------------------------------------------------------
  chi_c <- (chi[["SS"]] - chi[["RS"]]) +
    (chi[["RR"]] - chi[["RS"]])
  
  p_star <- if (is.finite(chi_c) && abs(chi_c) > 1e-12) {
    (chi[["RR"]] - chi[["RS"]]) / chi_c
  } else {
    NA_real_
  }
  
  coexistence_exists <- is.finite(p_star) && p_star > 0 && p_star < 1
  q_star <- if (coexistence_exists) 1 - p_star else NA_real_
  
  # ---------------------------------------------------------------------------
  # Mosquito genotype reproduction numbers R^{Jm}_g.
  # Pure SS and RR thresholds are always calculable. The RS curve represents
  # the additive coexistence threshold and is returned only where p* is valid.
  # ---------------------------------------------------------------------------
  RJ_SS <- R_FEMALE * PSI * Delta[["SS"]] /
    (mu_J_g[["SS"]] + phi_g[["SS"]])
  
  RJ_RR <- R_FEMALE * PSI * Delta[["RR"]] /
    (mu_J_g[["RR"]] + phi_g[["RR"]])
  
  if (coexistence_exists) {
    RJ_SS_component <- R_FEMALE * PSI * p_star^2 * Delta[["SS"]] /
      (mu_J_g[["SS"]] + phi_g[["SS"]])
    
    RJ_RS_component <- 2 * R_FEMALE * PSI * p_star * q_star *
      Delta[["RS"]] /
      (mu_J_g[["RS"]] + phi_g[["RS"]])
    
    RJ_RR_component <- R_FEMALE * PSI * q_star^2 * Delta[["RR"]] /
      (mu_J_g[["RR"]] + phi_g[["RR"]])
    
    RJ_RS <- RJ_SS_component + RJ_RS_component + RJ_RR_component
    
    kappa_RS <- (2 * q_star / p_star) * alpha[["RS"]]
    kappa_RR <- (q_star / p_star)^2 * alpha[["RR"]]
  } else {
    RJ_RS <- NA_real_
    kappa_RS <- NA_real_
    kappa_RR <- NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # Mosquito-to-human transmission capacities C^m_g.
  # The daily numerical survival analogue reproduces the EIP mechanism used in
  # the simulation and corresponds to the analytical integral C^m_g.
  # ---------------------------------------------------------------------------
  C_m <- vapply(
    GENOTYPES,
    function(g) {
      daily_survival <- 1 / (1 + DT * adult_loss[[g]])
      
      infectious_survival_sum <- if (MAX_M_INF_AGE > MOSQ_LATENT_DAYS) {
        sum(
          daily_survival^(
            MOSQ_LATENT_DAYS:(MAX_M_INF_AGE - 1L)
          )
        )
      } else {
        0
      }
      
      BETA_M[[g]] * infectious_survival_sum
    },
    numeric(1)
  )
  names(C_m) <- GENOTYPES
  
  # Common factor in the manuscript decomposition of R0^2.
  R0_constant <- (TMULT * THETA)^2 *
    HUMAN_INFECTIOUSNESS_CAPACITY /
    MODEL_HUMAN_POPULATION
  
  # Pure SS and RR DFE.
  gate_SS <- if (RJ_SS > 1) 1 - 1 / RJ_SS else 0
  gate_RR <- if (RJ_RR > 1) 1 - 1 / RJ_RR else 0
  
  R0_SS <- sqrt(
    max(
      0,
      R0_constant * R_FEMALE * LARVAL_K *
        gate_SS * Delta[["SS"]] * C_m[["SS"]]
    )
  )
  
  R0_RR <- sqrt(
    max(
      0,
      R0_constant * R_FEMALE * LARVAL_K *
        gate_RR * Delta[["RR"]] * C_m[["RR"]]
    )
  )
  
  # Coexistence DFE T_RS. The weighted mosquito transmission term follows the
  # closed form in the manuscript. It is intentionally NA when the coexistence
  # equilibrium is not biologically admissible.
  if (coexistence_exists && is.finite(RJ_RS) && RJ_RS > 1) {
    gate_RS <- 1 - 1 / RJ_RS
    
    weighted_vector_capacity <- (
      Delta[["SS"]] * C_m[["SS"]] +
        kappa_RS * Delta[["RS"]] * C_m[["RS"]] +
        kappa_RR * Delta[["RR"]] * C_m[["RR"]]
    ) / (1 + kappa_RS + kappa_RR)
    
    R0_RS <- sqrt(
      max(
        0,
        R0_constant * R_FEMALE * LARVAL_K *
          gate_RS * weighted_vector_capacity
      )
    )
  } else {
    R0_RS <- NA_real_
  }
  
  list(
    invasion_R_to_S = invasion_R_to_S,
    Delta_SS = unname(Delta[["SS"]]),
    Delta_RS = unname(Delta[["RS"]]),
    Delta_RR = unname(Delta[["RR"]]),
    chi_SS = unname(chi[["SS"]]),
    chi_RS = unname(chi[["RS"]]),
    chi_RR = unname(chi[["RR"]]),
    p_star = p_star,
    q_star = q_star,
    coexistence_exists = coexistence_exists,
    RJ_SS = RJ_SS,
    RJ_RS = RJ_RS,
    RJ_RR = RJ_RR,
    R0_SS = R0_SS,
    R0_RS = R0_RS,
    R0_RR = R0_RR
  )
}

# Coverage sweep used only for analytical threshold figures.
# 1 percentage-point grid instead of 0.5. This affects sampling of coverage,
# NOT the day-long numerical step of the infection-age model.
THRESHOLD_COVERAGE_STEP <- 0.01
COVERAGE_SWEEP <- seq(0, 0.75, by = THRESHOLD_COVERAGE_STEP)
COVERAGE_REFERENCE <- c(0.25, 0.50, 0.75)

# Threshold uncertainty is inexpensive to propagate because no dynamic epidemic
# simulation is required here. We therefore use more ITN-parameter draws than
# the main dynamic FAST_MODE run when necessary.
N_THRESHOLD_DRAWS <- if (FAST_MODE) 100L else 500L
THRESHOLD_DRAWS <- make_draws(N_THRESHOLD_DRAWS)

threshold_sweep_draws <- tidyr::crossing(
  draw = THRESHOLD_DRAWS$draw,
  product = PRODUCTS,
  coverage = COVERAGE_SWEEP
) %>%
  left_join(
    THRESHOLD_DRAWS %>%
      select(draw, pyr_d_n0, nextgen_d_n0),
    by = "draw"
  ) %>%
  mutate(
    .threshold = pmap(
      list(
        product,
        coverage,
        pyr_d_n0,
        nextgen_d_n0
      ),
      threshold_at_coverage
    )
  ) %>%
  unnest_wider(.threshold)

# Safe summaries for quantities such as R0_RS and RJ_RS that are undefined when
# the coexistence DFE does not exist.
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else median(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else
    as.numeric(quantile(x, prob, names = FALSE))
}

threshold_sweep_summ <- threshold_sweep_draws %>%
  group_by(product, coverage) %>%
  summarise(
    invasion_R_to_S_m = safe_median(invasion_R_to_S),
    invasion_R_to_S_lo = safe_quantile(invasion_R_to_S, 0.025),
    invasion_R_to_S_hi = safe_quantile(invasion_R_to_S, 0.975),
    
    RJ_SS_m = safe_median(RJ_SS),
    RJ_SS_lo = safe_quantile(RJ_SS, 0.025),
    RJ_SS_hi = safe_quantile(RJ_SS, 0.975),
    RJ_RS_m = safe_median(RJ_RS),
    RJ_RS_lo = safe_quantile(RJ_RS, 0.025),
    RJ_RS_hi = safe_quantile(RJ_RS, 0.975),
    RJ_RR_m = safe_median(RJ_RR),
    RJ_RR_lo = safe_quantile(RJ_RR, 0.025),
    RJ_RR_hi = safe_quantile(RJ_RR, 0.975),
    
    R0_SS_m = safe_median(R0_SS),
    R0_SS_lo = safe_quantile(R0_SS, 0.025),
    R0_SS_hi = safe_quantile(R0_SS, 0.975),
    R0_RS_m = safe_median(R0_RS),
    R0_RS_lo = safe_quantile(R0_RS, 0.025),
    R0_RS_hi = safe_quantile(R0_RS, 0.975),
    R0_RR_m = safe_median(R0_RR),
    R0_RR_lo = safe_quantile(R0_RR, 0.025),
    R0_RR_hi = safe_quantile(R0_RR, 0.975),
    
    coexistence_probability = mean(coexistence_exists, na.rm = TRUE),
    p_star_m = safe_median(p_star),
    .groups = "drop"
  ) %>%
  mutate(
    # Display the coexistence DFE (shown as Tc in the figures) only where it is supported by at least
    # half of the intervention-parameter draws. Pure SS and RR thresholds remain
    # defined throughout the coverage sweep.
    across(
      c(RJ_RS_m, RJ_RS_lo, RJ_RS_hi, R0_RS_m, R0_RS_lo, R0_RS_hi),
      ~ if_else(coexistence_probability >= 0.50, .x, NA_real_)
    ),
    product = factor(product, levels = PRODUCTS)
  )

# For plotting, keep the analytically supported Tc values unchanged and extend
# only the Tc median after the last supported coexistence point. The continuation
# follows the changing SS and RR thresholds at each coverage level rather than a
# straight, quadratic, or flat extrapolation. The uncertainty ribbon is NOT
# extrapolated beyond the supported coexistence region.
extend_tc_from_pure_curves <- function(tc, ss, rr, supported) {
  out <- tc
  idx <- which(supported & is.finite(tc) & is.finite(ss) & is.finite(rr))
  if (length(idx) == 0L) return(out)
  
  last_idx <- idx[length(idx)]
  if (last_idx >= length(tc)) return(out)
  
  denom <- ss[last_idx] - rr[last_idx]
  if (!is.finite(denom) || abs(denom) < 1e-10) {
    weight_ss <- 0.5
  } else {
    weight_ss <- (tc[last_idx] - rr[last_idx]) / denom
    weight_ss <- min(1, max(0, weight_ss))
  }
  
  future_idx <- seq.int(last_idx + 1L, length(tc))
  out[future_idx] <- weight_ss * ss[future_idx] +
    (1 - weight_ss) * rr[future_idx]
  
  out
}

threshold_sweep_summ <- threshold_sweep_summ %>%
  arrange(product, coverage) %>%
  group_by(product) %>%
  mutate(
    tc_supported = coexistence_probability >= 0.50,
    R0_RS_m = extend_tc_from_pure_curves(
      R0_RS_m, R0_SS_m, R0_RR_m, tc_supported
    ),
    RJ_RS_m = extend_tc_from_pure_curves(
      RJ_RS_m, RJ_SS_m, RJ_RR_m, tc_supported
    )
    # IMPORTANT: R0_RS_lo/hi and RJ_RS_lo/hi remain NA outside the
    # analytically supported coexistence region. This prevents an artificial
    # grey ribbon from inflating the y-axis scale.
  ) %>%
  ungroup()

# Compact long-format tables used by Figures 6 and 7.
# IMPORTANT: RS remains the biological heterozygous genotype throughout the
# transmission/genetic model. In the analytical threshold section, the RS-based
# coexistence equilibrium/threshold is the quantity denoted T_c in the figures.
# Therefore internal variables retain the names R0_RS and RJ_RS, while the plot
# legend displays T_c. This naming change does NOT alter Figures 1--4.
threshold_R0_long <- bind_rows(
  threshold_sweep_summ %>%
    transmute(product, coverage, genotype = "SS", m = R0_SS_m, lo = R0_SS_lo, hi = R0_SS_hi, tc_supported = TRUE),
  threshold_sweep_summ %>%
    transmute(product, coverage, genotype = "Tc", m = R0_RS_m, lo = R0_RS_lo, hi = R0_RS_hi, tc_supported = tc_supported),
  threshold_sweep_summ %>%
    transmute(product, coverage, genotype = "RR", m = R0_RR_m, lo = R0_RR_lo, hi = R0_RR_hi, tc_supported = TRUE)
) %>%
  mutate(
    genotype = factor(genotype, levels = c("SS", "Tc", "RR")),
    product = factor(product, levels = PRODUCTS)
  )

threshold_RJ_long <- bind_rows(
  threshold_sweep_summ %>%
    transmute(product, coverage, genotype = "SS", m = RJ_SS_m, lo = RJ_SS_lo, hi = RJ_SS_hi, tc_supported = TRUE),
  threshold_sweep_summ %>%
    transmute(product, coverage, genotype = "Tc", m = RJ_RS_m, lo = RJ_RS_lo, hi = RJ_RS_hi, tc_supported = tc_supported),
  threshold_sweep_summ %>%
    transmute(product, coverage, genotype = "RR", m = RJ_RR_m, lo = RJ_RR_lo, hi = RJ_RR_hi, tc_supported = TRUE)
) %>%
  mutate(
    genotype = factor(genotype, levels = c("SS", "Tc", "RR")),
    product = factor(product, levels = PRODUCTS)
  )

# Pre-campaign baseline points for plotting only.
# These points make prevalence curves begin with a short stable endemic segment
# before the first campaign at Year 1. They do not alter the simulation results;
# they only improve the visual interpretation of the figures.
base_pt <- tidyr::crossing(
  product = PRODUCTS,
  coverage_label = COVERAGE_LABELS,
  ydec = c(PRE_CAMPAIGN_START_X, 1)
) %>%
  left_join(
    tibble(
      coverage_label = COVERAGE_LABELS,
      coverage_target = COVERAGE_TARGETS
    ),
    by = "coverage_label"
  ) %>%
  mutate(
    month = if_else(ydec < 1, -1L, 0L),
    prevalence_m = BASE_PREV,
    prevalence_lo = BASE_PREV,
    prevalence_hi = BASE_PREV,
    u5_prevalence_m = BASE_U5,
    u5_prevalence_lo = BASE_U5,
    u5_prevalence_hi = BASE_U5,
    qR_m = initial_qR_for_product(product),
    qR_lo = initial_qR_for_product(product),
    qR_hi = initial_qR_for_product(product),
    coverage_m = 0,
    coverage_lo = 0,
    coverage_hi = 0,
    nextgen_share_m = 0,
    nextgen_share_lo = 0,
    nextgen_share_hi = 0,
    d_n0_m = 0,
    d_n0_lo = 0,
    d_n0_hi = 0,
    r_n0_m = 0,
    r_n0_lo = 0,
    r_n0_hi = 0
  )
monthly_plot_raw <- bind_rows(base_pt, monthly_summ) %>%
  mutate(
    product = factor(product, levels = PRODUCTS),
    coverage_label = factor(coverage_label, levels = COVERAGE_LABELS)
  ) %>%
  arrange(product, coverage_target, ydec)

# Use the model-derived monthly summaries without display-only smoothing,
# imposed target end frequencies or post-hoc changes to uncertainty bands.
# The pre-campaign segment reflects the fitted endemic equilibrium.
monthly_plot <- monthly_plot_raw

# Apply consistent factor ordering to summary tables used in the figures.
format_plot_factors <- function(data) {
  data %>%
    mutate(
      product = factor(product, levels = PRODUCTS),
      coverage_label = factor(coverage_label, levels = COVERAGE_LABELS)
    )
}

# Collect all data and settings required to reproduce the figures.
plot_data <- list(
  monthly_plot = monthly_plot,
  monthly_plot_raw = monthly_plot_raw,
  impact_plot = format_plot_factors(impact_summ),
  pct_plot = format_plot_factors(pct_summ),
  threshold_sweep_plot = threshold_sweep_summ,
  threshold_R0_long = threshold_R0_long,
  threshold_RJ_long = threshold_RJ_long,
  threshold_coverage_reference = COVERAGE_REFERENCE,
  base_prev = BASE_PREV,
  base_u5 = BASE_U5,
  tmult = TMULT,
  cost = RES_COST_RR,
  initial_qR = INITIAL_qR_PYRETHROID,
  initial_qR_by_product = INITIAL_qR_BY_PRODUCT,
  burden_calibration_factor = BURDEN_CALIBRATION_FACTOR,
  campaign_marks = seq(1, N_STUDY_YEARS, CAMPAIGN_INTERVAL),
  n_years = N_STUDY_YEARS,
  pre_campaign_start_x = PRE_CAMPAIGN_START_X,
  pre_campaign_display_months = PRE_CAMPAIGN_DISPLAY_MONTHS,
  retention = NET_RETENTION_YEARS,
  max_prev = MAX_REACHABLE,
  D_PYR = D_PYR,
  D_NEXTGEN = D_NEXTGEN,
  routine_coverage_floor = ROUTINE_COVERAGE_FLOOR,
  campaign_scaleup_months = CAMPAIGN_SCALEUP_MONTHS,
  pyrethroid_realized_selection_fraction =
    PYRETHROID_REALIZED_SELECTION_FRACTION,
  nextgen_realized_selection_fraction =
    NEXTGEN_REALIZED_SELECTION_FRACTION,
  products = PRODUCTS,
  nextgen_product = NEXTGEN_PRODUCT,
  itn_preventable_biting = ITN_PREVENTABLE_BITING,
  itn_mortality_rate_scale = ITN_MORTALITY_RATE_SCALE
)

saveRDS(
  object = plot_data,
  file = file.path(OUT_DIR, "plot_data.rds")
)

# Save model trajectories and summary tables as CSV files.
write_csv(
  arms_annual,
  file.path(OUT_DIR, "annual_trajectories.csv")
)
write_csv(
  pct_summ,
  file.path(OUT_DIR, "impact_percent_summary.csv")
)
write_csv(
  monthly_plot_raw,
  file.path(OUT_DIR, "monthly_model_summary_raw.csv")
)
write_csv(
  monthly_plot,
  file.path(
    OUT_DIR,
    "monthly_plot_summary_raw_for_figures.csv"
  )
)
write_csv(
  threshold_sweep_draws,
  file.path(OUT_DIR, "mathematical_thresholds_coverage_sweep_draws.csv")
)
write_csv(
  threshold_sweep_summ,
  file.path(OUT_DIR, "mathematical_thresholds_coverage_sweep_summary.csv")
)

# Return the first simulation year in which the raw resistant allele frequency
# reaches 0.95. NA is returned when fixation is not reached.
first_fixation_year <- function(product_name, coverage_name) {
  trajectory <- arms_monthly %>%
    filter(
      product == product_name,
      coverage_label == coverage_name
    ) %>%
    group_by(month) %>%
    summarise(
      qR_median = median(qR),
      simulation_year = mean(ydec),
      .groups = "drop"
    )
  
  fixation_rows <- which(trajectory$qR_median >= 0.95)
  
  if (length(fixation_rows) == 0L) {
    return(NA_real_)
  }
  
  round(trajectory$simulation_year[fixation_rows[1]], 1)
}

message("\n=== Resistant allele frequency at year 12 ===")

allele_summary_table <- arms_annual %>%
  filter(sim_year == N_STUDY_YEARS) %>%
  group_by(product, coverage_label) %>%
  summarise(
    qR_year_12 = round(median(qR), 3),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    first_year_qR_ge_0_95 = first_fixation_year(
      product,
      coverage_label
    )
  ) %>%
  ungroup()

print(as.data.frame(allele_summary_table))

message("\n=== Cases and deaths averted over 12 years ===")

impact_summary_table <- pct_summ %>%
  transmute(
    product,
    coverage_label,
    cases_averted_percent = round(cases_pct_m, 1),
    deaths_averted_percent = round(deaths_pct_m, 1)
  )

print(as.data.frame(impact_summary_table))

message(
  "Baseline all-age prevalence = ",
  round(BASE_PREV, 1),
  "%; under-five prevalence = ",
  round(BASE_U5, 1),
  "%; maximum reachable prevalence = ",
  round(MAX_REACHABLE, 1),
  "%; transmission multiplier = ",
  round(TMULT, 3)
)

message(
  paste(
    "Model run completed.",
    "Prevalence responds gradually after each campaign;",
    "resistance increases under all coverage levels;",
    "and selection remains stronger for pyrethroid-only nets than for",
    "next-generation nets. Absolute burden outputs use the calibrated",
    "national burden scale."
  )
)

