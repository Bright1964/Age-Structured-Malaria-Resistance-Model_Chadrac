# Figure 2: single retained daily mechanistic provincial model.
# This script should run after the main and supplement scripts.

# =============================================================================
# DRC: RESISTANT-ALLELE FREQUENCY AND PROVINCIAL ITN CAMPAIGNS
# ONE FIGURE ONLY: A (9 provinces), B (9 provinces), C (8 provinces) | FAST
# =============================================================================
# Campaigns come from Excel: a value of 1 indicates a campaign in that year.
# Coverage target = 50% for every campaign; retention follows the supplied model.
# q_R is a raw mechanistic prediction (NOT the observed provincial frequency).
# All provinces use the same initial q_R and biological parameters; only the
# campaign timing and resulting coverage differ between provinces.
#
# IMPORTANT: Excel does not contain campaign product types. The use of a
# 2019 PBO cutoff below is a WORKING ASSUMPTION, NOT an observed provincial
# distribution history. Set PBO_FIRST_YEAR <- NA_integer_ to disable it, or
# replace the classification with verified province-specific product data.
# The original model's "Next-generation nets" parameters are used as a proxy
# for PBO nets when the assumption is enabled; this is not a calibrated PBO model.
# =============================================================================

RUN_START_TIME <- proc.time()[["elapsed"]]

packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "scales", "readr")
missing_pkgs <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Install missing R packages before running: ",
       paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ---------------------- 1. USER SETTINGS -------------------------------------
dir.create("cache", recursive = TRUE, showWarnings = FALSE)
dir.create("output", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
EXCEL_FILE <- file.path("data", "bednet_distribution_by_province_2008_2022.xlsx")
OUTPUT_PDF <- file.path("figures", "Figure_2_provincial_resistance.pdf")

# Speed options: keeps the same daily mechanistic model and calibration settings.
# On Windows, PSOCK parallel workers avoid the unsupported multicore fork mode.
USE_PARALLEL <- TRUE
CPU_AVAILABLE <- parallel::detectCores(logical = FALSE)
if (!is.finite(CPU_AVAILABLE) || CPU_AVAILABLE < 1L) CPU_AVAILABLE <- 2L
N_WORKERS <- max(1L, min(4L, as.integer(CPU_AVAILABLE) - 1L))

# Saving numerical results makes repeat runs (e.g. figure edits) very fast.
USE_CACHE <- TRUE
FORCE_RECALIBRATE <- FALSE
FORCE_PROVINCE_SIMULATIONS <- FALSE
SPINUP_CACHE_FILE <- file.path("cache", "resistance_model_spinup_cache.rds")
PROVINCE_CACHE_FILE <- file.path("cache", "resistance_province_trajectories_cache.rds")

STUDY_START_YEAR <- 2009L
STUDY_END_YEAR   <- 2024L
STUDY_YEARS <- STUDY_START_YEAR:STUDY_END_YEAR
N_STUDY_YEARS <- length(STUDY_YEARS)
DAYS_PER_YEAR <- 365L
# Daily step is retained; the infection-age array advances exactly one day.
DT <- 1.0
SIM_LENGTH <- N_STUDY_YEARS * DAYS_PER_YEAR
N_MONTHS <- 12L * N_STUDY_YEARS

# Target coverage is exactly 50% at the campaign peak (not 50% every day).
COVERAGE_TARGET <- 0.50
PEAK_MULTIPLIER <- 1.00

# 2019 is a modelling/graphics assumption, NOT a fact in the Excel sheet.
# Use NA_integer_ if you want to mark all recorded campaigns as pyrethroid.
PBO_FIRST_YEAR <- 2019L

# Preserve the numerical settings of the supplied mechanistic model.
PUBLICATION_COMMON_INITIAL_qR <- 0.05
TARGET_BASELINE_PREVALENCE <- 0.50
PLOT_BASELINE_PREVALENCE <- 50
MODEL_HUMAN_POPULATION <- 1e5
NATIONAL_POP <- 70e6
SPINUP_YEARS <- 20L
CAL_SPINUP_YEARS <- 15L
CAL_ITERS <- 10L

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

# Compile the identical daily model into R bytecode. No C++/Rtools needed.
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

# Cache is reused ONLY when the mathematical input parameters agree.
# Change MODEL_CACHE_VERSION whenever the equations in this script are revised.
MODEL_CACHE_VERSION <- 2L
MODEL_SIGNATURE <- list(
  version = MODEL_CACHE_VERSION,
  target = TARGET_BASELINE_PREVALENCE,
  spinup_years = SPINUP_YEARS,
  cal_spinup_years = CAL_SPINUP_YEARS,
  cal_iters = CAL_ITERS,
  age_bins = AGE_BINS_YEARS,
  recovery = RECOVERY_DAYS,
  cfr = CFR_AGE,
  p_clin = P_CLIN_AGE,
  mu_h = MU_H,
  kh = K_H,
  theta = THETA,
  beta_m = BETA_M,
  phi_m = PHI_M,
  mu_jm = MU_JM,
  mu_am = MU_AM,
  psi = PSI,
  larval_k = LARVAL_K,
  eir = NATIONAL_EIR,
  genotype_par = c(INITIAL_qR, RES_h, RES_u,
                   RES_COST_RR, RES_COST_RS_FRACTION),
  effects = c(ITN_PREVENTABLE_BITING,
              ITN_MORTALITY_RATE_SCALE,
              PYRETHROID_REALIZED_SELECTION_FRACTION,
              NEXTGEN_REALIZED_SELECTION_FRACTION),
  H2M = H2M,
  population = MODEL_HUMAN_POPULATION
)

spinup_cache_ok <- FALSE
if (USE_CACHE && !FORCE_RECALIBRATE &&
    file.exists(SPINUP_CACHE_FILE)) {
  prior_spinup <- tryCatch(readRDS(SPINUP_CACHE_FILE),
                          error = function(e) NULL)
  spinup_cache_ok <- is.list(prior_spinup) &&
    identical(prior_spinup$signature, MODEL_SIGNATURE) &&
    is.list(prior_spinup$SPIN) &&
    is.numeric(prior_spinup$TMULT) &&
    length(prior_spinup$TMULT) == 1L &&
    is.finite(prior_spinup$TMULT)
}

if (spinup_cache_ok) {
  message("Loading validated spin-up/calibration cache (no recalculation).")
  SPIN <- prior_spinup$SPIN
  TMULT <- prior_spinup$TMULT
  TARGET_USED <- prior_spinup$TARGET_USED
  MAX_REACHABLE <- prior_spinup$MAX_REACHABLE
} else {
  message("Calibrating initial model (first run takes longer)...")
  MAX_REACHABLE <- endemic_prev(1e3) * 100
  TARGET_USED <- min(TARGET_BASELINE_PREVALENCE,
                     0.95 * MAX_REACHABLE / 100)
  TMULT <- calibrate_tmult(TARGET_USED)
  SPIN <- run_spinup(TMULT, SPINUP_YEARS)
  if (USE_CACHE) {
    saveRDS(list(signature = MODEL_SIGNATURE, SPIN = SPIN,
                 TMULT = TMULT, TARGET_USED = TARGET_USED,
                 MAX_REACHABLE = MAX_REACHABLE),
            file = SPINUP_CACHE_FILE, compress = FALSE)
  }
}

message("Transmission multiplier = ", round(TMULT, 4),
        "; baseline all-age prevalence = ",
        round(SPIN$prevalence, 2), "%")

# =============================================================================
# 2. PROVINCIAL CAMPAIGNS FROM EXCEL (ACTUAL 0/1 SCHEDULES)
# =============================================================================
if (!file.exists(EXCEL_FILE)) {
  # Also search alongside the current sourced .R file, not only getwd().
  source_files <- unlist(lapply(sys.frames(), function(f) {
    if (is.null(f$ofile)) character(0) else as.character(f$ofile)
  }), use.names = FALSE)
  if (length(source_files)) {
    candidate_paths <- file.path(
      dirname(normalizePath(source_files, mustWork = FALSE)),
      basename(EXCEL_FILE)
    )
    existing <- candidate_paths[file.exists(candidate_paths)]
    if (length(existing)) EXCEL_FILE <- existing[1L]
  }
}
if (!file.exists(EXCEL_FILE)) {
  stop("Excel file not found: ", EXCEL_FILE,
       ". Place the .xlsx next to the .R file or set EXCEL_FILE to its full path.")
}

campaign_wide <- readxl::read_excel(EXCEL_FILE)
province_col <- names(campaign_wide)[1]
if (tolower(province_col) != "province") {
  stop("The first Excel column must be 'province'.")
}
if (nrow(campaign_wide) != 26L) {
  stop("Expected 26 provinces; found ", nrow(campaign_wide))
}
if (anyDuplicated(campaign_wide$province)) {
  stop("Duplicate province names were found in Excel.")
}

required_year_cols <- as.character(STUDY_YEARS)
missing_years <- setdiff(required_year_cols, names(campaign_wide))
if (length(missing_years) > 0L) {
  stop("Missing years in Excel: ", paste(missing_years, collapse = ", "))
}

# Keep EXACTLY the order of provinces in Excel (do not sort alphabetically).
province_info <- tibble(
  province = as.character(campaign_wide$province),
  panel = rep(c("A", "B", "C"), times = c(9L, 9L, 8L))
) %>%
  mutate(
    province_key = paste(panel, province, sep = "::")
  )

campaign_long <- campaign_wide %>%
  select(province, all_of(required_year_cols)) %>%
  mutate(province = as.character(province)) %>%
  pivot_longer(
    cols = all_of(required_year_cols),
    names_to = "year",
    values_to = "distributed"
  ) %>%
  mutate(
    year = as.integer(year),
    distributed = as.numeric(distributed)
  )

if (anyNA(campaign_long$distributed) ||
    !all(campaign_long$distributed %in% c(0, 1))) {
  stop("Campaign values must be 0 or 1, with no missing entries.")
}

campaign_long <- campaign_long %>%
  filter(distributed == 1) %>%
  left_join(province_info, by = "province") %>%
  mutate(
    campaign_type = if (is.na(PBO_FIRST_YEAR)) {
      "Pyrethroid"
    } else {
      if_else(year >= PBO_FIRST_YEAR, "PBO (assumed)", "Pyrethroid")
    },
    province_key = factor(
      province_key,
      levels = rev(province_info$province_key)
    )
  ) %>%
  arrange(match(province, province_info$province), year)

# =============================================================================
# 3. COVERAGE C_B(t): PROVINCE-SPECIFIC CAMPAIGNS, 50% PEAK
# =============================================================================
# Adaptation of build_coverage() from the supplied R code; only the campaign
# schedule is made province-specific. The six-month build-up, waning and
# 10% routine floor match the original model structure.
build_province_coverage <- function(campaign_years, target = COVERAGE_TARGET) {
  campaign_days <- as.integer(
    (sort(unique(campaign_years)) - STUDY_START_YEAR) * DAYS_PER_YEAR + 1L
  )
  campaign_days <- campaign_days[
    campaign_days >= 1L & campaign_days <= SIM_LENGTH
  ]

  cov <- numeric(SIM_LENGTH)
  if (length(campaign_days) == 0L) return(cov)

  is_campaign_day <- rep(FALSE, SIM_LENGTH)
  is_campaign_day[campaign_days] <- TRUE

  peak <- min(0.95, max(ROUTINE_COVERAGE_FLOOR,
                        target * PEAK_MULTIPLIER))
  routine_floor <- min(ROUTINE_COVERAGE_FLOOR, 0.80 * target)
  decay_rate <- log(2) / (NET_RETENTION_YEARS * DAYS_PER_YEAR)

  last_campaign_day <- NA_integer_
  campaign_start_coverage <- 0

  for (day in seq_len(SIM_LENGTH)) {
    if (is_campaign_day[day]) {
      last_campaign_day <- day
      campaign_start_coverage <- if (day == 1L) 0 else cov[day - 1L]
    }

    if (is.na(last_campaign_day)) next

    days_since_campaign <- day - last_campaign_day
    if (days_since_campaign <= CAMPAIGN_SCALEUP_DAYS) {
      fraction <- days_since_campaign / max(1, CAMPAIGN_SCALEUP_DAYS)
      smoothstep <- 3 * fraction^2 - 2 * fraction^3
      cov[day] <- campaign_start_coverage +
        (peak - campaign_start_coverage) * smoothstep
    } else {
      days_since_peak <- days_since_campaign - CAMPAIGN_SCALEUP_DAYS
      cov[day] <- routine_floor +
        (peak - routine_floor) * exp(-decay_rate * days_since_peak)
    }
  }

  pmin(0.95, pmax(0, cov))
}

# Product carried forward until the next campaign, for the modelled net effect.
# PBO is proxied using the "next-generation" coefficients of the original code.
product_by_day <- function(campaign_years, campaign_types) {
  type <- rep("Pyrethroid", SIM_LENGTH)
  if (length(campaign_years) == 0L) return(type)

  campaign_days <- as.integer(
    (campaign_years - STUDY_START_YEAR) * DAYS_PER_YEAR + 1L
  )
  for (i in seq_along(campaign_days)) {
    end_day <- if (i < length(campaign_days)) {
      campaign_days[i + 1L] - 1L
    } else {
      SIM_LENGTH
    }
    type[campaign_days[i]:end_day] <- campaign_types[i]
  }
  type
}

# =============================================================================
# 4. ONE PROVINCE = ONE MECHANISTIC SIMULATION
# =============================================================================
# Unlike display-only postprocessing, the allele frequencies here are read
# directly from the female mosquito genotypes in step_model().
month_end_days <- as.integer(round(
  seq_len(N_MONTHS) * DAYS_PER_YEAR / 12
))

# Use base R in workers: no package loading, minimal communication overhead.
simulate_task <- function(task) {
  coverage <- build_province_coverage(task$years)
  products <- product_by_day(task$years, task$types)
  state <- reset_state_allele_frequency(
    SPIN, PUBLICATION_COMMON_INITIAL_qR
  )

  q_monthly <- numeric(N_MONTHS)
  next_month <- 1L

  for (day in seq_len(SIM_LENGTH)) {
    is_pbo <- identical(products[day], "PBO (assumed)")

    state <- step_model(
      state = state,
      coverage = coverage[day],
      d_n0 = if (is_pbo) 0.58 else 0.48,
      r_n0 = if (is_pbo) 0.52 else 0.48,
      d_g = if (is_pbo) D_NEXTGEN else D_PYR,
      selection_fraction = if (is_pbo) {
        NEXTGEN_REALIZED_SELECTION_FRACTION
      } else {
        PYRETHROID_REALIZED_SELECTION_FRACTION
      },
      day = day,
      include_cost = TRUE,
      tmult = TMULT
    )

    if (next_month <= N_MONTHS && day == month_end_days[next_month]) {
      q_monthly[next_month] <- state$qR
      next_month <- next_month + 1L
    }
  }

  if (next_month != N_MONTHS + 1L) {
    stop("Not all monthly states were recorded for ", task$province)
  }

  data.frame(
    province = task$province,
    month = seq_len(N_MONTHS),
    year_decimal = STUDY_START_YEAR + (seq_len(N_MONTHS) - 0.5) / 12,
    qR = q_monthly,
    stringsAsFactors = FALSE
  )
}
simulate_task <- compiler::cmpfun(simulate_task)

# All province schedules are taken directly from the Excel (never from a
# national 3-year cycle). Repeated schedules are simulated only once.
tasks <- lapply(province_info$province, function(prov) {
  events <- campaign_long[campaign_long$province == prov,
                          c("year", "campaign_type"), drop = FALSE]
  events <- events[order(events$year), , drop = FALSE]
  list(province = prov, years = events$year, types = events$campaign_type)
})

# Include campaign types in the key, not just campaign years.
schedule_keys <- vapply(tasks, function(task) {
  paste(paste(task$years, task$types, sep = ":"), collapse = ";")
}, character(1))
unique_indices <- which(!duplicated(schedule_keys))
unique_tasks <- tasks[unique_indices]
unique_keys <- schedule_keys[unique_indices]

PROVINCE_SIGNATURE <- list(
  model = MODEL_SIGNATURE,
  campaigns = schedule_keys,
  coverage_target = COVERAGE_TARGET,
  peak_multiplier = PEAK_MULTIPLIER,
  retention = NET_RETENTION_YEARS,
  routine_floor = ROUTINE_COVERAGE_FLOOR,
  scaleup_days = CAMPAIGN_SCALEUP_DAYS,
  pbo_cutoff = PBO_FIRST_YEAR,
  nextgen_restore = NEXTGEN_RESTORE,
  time = c(STUDY_START_YEAR, STUDY_END_YEAR),
  q_initial = PUBLICATION_COMMON_INITIAL_qR,
  spinup = c(TMULT, SPIN$prevalence)
)

province_cache_ok <- FALSE
if (USE_CACHE && !FORCE_PROVINCE_SIMULATIONS &&
    file.exists(PROVINCE_CACHE_FILE)) {
  prior_provinces <- tryCatch(readRDS(PROVINCE_CACHE_FILE),
                              error = function(e) NULL)
  province_cache_ok <- is.list(prior_provinces) &&
    identical(prior_provinces$signature, PROVINCE_SIGNATURE) &&
    is.data.frame(prior_provinces$monthly_resistance) &&
    nrow(prior_provinces$monthly_resistance) == 26L * N_MONTHS
}

if (province_cache_ok) {
  message("Loading validated 26-province trajectory cache (no simulations).")
  monthly_resistance <- prior_provinces$monthly_resistance
} else {
  message("Simulating ", length(unique_tasks),
          " distinct schedules for 26 provinces using ", N_WORKERS,
          " worker(s)...")

  # Required by the model functions, their helper functions and the solver.
  model_worker_vars <- c(
    "SPIN", "TMULT", "D_PYR", "D_NEXTGEN", "N_MONTHS",
    "SIM_LENGTH", "STUDY_START_YEAR", "DAYS_PER_YEAR", "DT",
    "COVERAGE_TARGET", "PEAK_MULTIPLIER", "NET_RETENTION_YEARS",
    "ROUTINE_COVERAGE_FLOOR", "CAMPAIGN_SCALEUP_DAYS",
    "PUBLICATION_COMMON_INITIAL_qR", "GENOTYPES", "N_G", "N_AGE",
    "MU_AM", "MU_JM", "MU_H", "GAMMA_H_AGE", "NU_H_AGE", "K_H",
    "MODEL_HUMAN_POPULATION", "AGEING_RATE", "MAX_H_INF_AGE",
    "MAX_M_INF_AGE", "MOSQ_LATENT_DAYS", "THETA", "BETA_M",
    "H2M", "PHI_M", "PSI", "R_FEMALE", "LARVAL_K", "UNDER5",
    "RES_COST_RR", "RES_COST_RS_FRACTION",
    "ITN_MORTALITY_RATE_SCALE", "ITN_PREVENTABLE_BITING",
    "NEXTGEN_REALIZED_SELECTION_FRACTION",
    "PYRETHROID_REALIZED_SELECTION_FRACTION", "month_end_days",
    "sanitize_positive", "normalise_vector", "sp", "norm_vec",
    "cost_vec", "step_model", "reset_state_allele_frequency",
    "build_province_coverage", "product_by_day", "simulate_task"
  )
  missing_vars <- setdiff(model_worker_vars, ls(envir = environment(),
                                                  all.names = TRUE))
  if (length(missing_vars) > 0L) {
    stop("Missing worker variables: ", paste(missing_vars, collapse = ", "))
  }

  nworkers <- min(N_WORKERS, length(unique_tasks))
  if (USE_PARALLEL && nworkers >= 2L) {
    cluster <- NULL
    unique_results <- tryCatch({
      cluster <- parallel::makeCluster(nworkers, type = "PSOCK")
      parallel::clusterExport(cluster, model_worker_vars,
                              envir = environment())
      parallel::parLapplyLB(cluster, unique_tasks, simulate_task)
    }, error = function(e) {
      message("Parallel computation unavailable: ", conditionMessage(e),
              ". Reverting to sequential simulations.")
      lapply(unique_tasks, simulate_task)
    }, finally = {
      if (!is.null(cluster)) parallel::stopCluster(cluster)
    })
  } else {
    unique_results <- lapply(unique_tasks, simulate_task)
  }

  names(unique_results) <- unique_keys
  monthly_resistance <- dplyr::bind_rows(lapply(seq_along(tasks), function(i) {
    result <- unique_results[[schedule_keys[i]]]
    result$province <- tasks[[i]]$province
    result
  }))

  if (USE_CACHE) {
    saveRDS(list(signature = PROVINCE_SIGNATURE,
                 monthly_resistance = monthly_resistance),
            file = PROVINCE_CACHE_FILE, compress = FALSE)
  }
}

monthly_resistance <- monthly_resistance %>%
  left_join(province_info, by = "province") %>%
  mutate(province_key = factor(
    province_key, levels = rev(province_info$province_key)
  ))

if (anyNA(monthly_resistance$qR) ||
    any(monthly_resistance$qR < 0 | monthly_resistance$qR > 1)) {
  stop("Invalid resistance allele frequency values in simulation results.")
}

# =============================================================================
# 5. EXACTLY ONE FIGURE WITH THREE PANELS: A (9), B (9), C (8)
# =============================================================================
# Use one common resistance colour scale [0,1] for all 26 provinces.
# The campaign symbols are located only in years that equal 1 in Excel.
figure_ABC <- ggplot(
  monthly_resistance,
  aes(x = year_decimal, y = province_key, fill = qR)
) +
  geom_tile(width = 1 / 12, height = 0.97) +
  geom_point(
    data = campaign_long,
    aes(
      x = year + 0.5,
      y = province_key,
      shape = campaign_type
    ),
    inherit.aes = FALSE,
    size = 1.8,
    colour = "black"
  ) +
  facet_wrap(
    ~ panel,
    nrow = 1,
    scales = "free_y",
    labeller = as_labeller(c(A = "A.", B = "B.", C = "C."))
  ) +
  scale_y_discrete(
    labels = function(x) sub("^[ABC]::", "", x)
  ) +
  scale_x_continuous(
    limits = c(STUDY_START_YEAR, STUDY_END_YEAR + 1),
    breaks = c(2009, 2012, 2015, 2018, 2021, 2024),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_gradientn(
    colours = c("#1b0c41", "#721f81", "#b63679", "#f1605d", "#feb078"),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = label_number(accuracy = 0.01),
    name = "Resistance\nallele frequency"
  ) +
  scale_shape_manual(
    values = c("Pyrethroid" = 17, "PBO (assumed)" = 16),
    breaks = c("Pyrethroid", "PBO (assumed)"),
    name = "Campaign"
  ) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 3)),
    fill = guide_colourbar(order = 2, barheight = grid::unit(3.5, "cm"))
  ) +
  labs(
    x = "Year",
    y = "Province",
    caption = if (is.na(PBO_FIRST_YEAR)) {
      "Campaign years: provincial Excel file. Raw mechanistic prediction; 50% target ITN coverage."
    } else {
      paste0(
        "Campaign years: provincial Excel file. PBO after ", PBO_FIRST_YEAR,
        " is an assumption (next-generation proxy); 50% target ITN coverage."
      )
    }
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15, hjust = 0),
    panel.border = element_rect(
      colour = "grey35", fill = NA, linewidth = 0.5
    ),
    panel.spacing.x = grid::unit(1.1, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.caption = element_text(size = 9, hjust = 0),
    plot.margin = margin(10, 10, 10, 10)
  )

# Export numerical source and campaign flags for reproducibility.
dir.create("output", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("cache", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(as.data.frame(monthly_resistance),
                 file.path("output", "Figure_2_monthly_provincial_resistance.csv"))
readr::write_csv(as.data.frame(campaign_long),
                 file.path("output", "Figure_2_campaign_markers.csv"))

# Exactly one publication figure (three panels A/B/C).
print(figure_ABC)

ggsave(
  filename = OUTPUT_PDF,
  plot = figure_ABC,
  device = "pdf",
  width = 20,
  height = 8.5,
  units = "in"
)

ggsave(
  filename = file.path("figures", "Figure_2_provincial_resistance.png"),
  plot = figure_ABC, width = 20, height = 8.5, units = "in", dpi = 300,
  bg = "white"
)
message("One three-panel figure saved: ", normalizePath(OUTPUT_PDF))
message("Total elapsed time: ",
        round(proc.time()[["elapsed"]] - RUN_START_TIME, 1), " seconds")


