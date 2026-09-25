# ==============================================================================
# 13-sequential-site.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Per Ellen Francis (September 2026):
#   - PRIMARY model now includes study site (fully adjusted + site + IPW)
#   - SECONDARY model without site (previous primary)
#   - Sequential covariate adjustment as a sensitivity analysis, to show
#     which covariates change the estimate between the minimally and fully
#     adjusted models, with particular attention to site
#   - Boston vs not Boston stratified models as a supplement
#
# Uses the corrected pipeline (07-10): 45,268 singleton pregnancies,
# 40,700 with an observed IH outcome, 20 imputations, Rubin's rules.
# Modified Poisson GEE, exchangeable correlation, clustered on MOMID.
#
# PART A. Sequential adjustment (ih1, unweighted so each step isolates
#         the covariate; IPW changed estimates by <= 0.01 in 10):
#   A1. Crude + each covariate individually (12 models per exposure)
#   A2. Cumulative blocks: crude → + demographics → + health/behaviour →
#       + reproductive/infant (= full, no site) → + site
#   Final rows: full + IPW (secondary) and full + site + IPW (primary)
# PART B. Boston (site 5) vs all other sites, fully adjusted + IPW
#         (other sites model also adjusts for site), ih1 and ih2
#
# Outputs (outputs/):
#   seq_individual.csv    A1, with % change in log RR vs crude
#   seq_cumulative.csv    A2 plus primary and secondary
#   boston_stratified.csv B
# No p-values.
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"   # folder with weighted_long_corrected.RDS
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE)
n_cores <- min(4, detectCores())   # set to 1 on Windows

# ------------------------------------------------------------------------------
# DATA (as 10)
# ------------------------------------------------------------------------------
weighted <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  mutate(
    smoking   = factor(smoking,   levels = c("Non-smoker","1-19/day","20+/day")),
    parity    = factor(parity,    levels = c("0","1","2","3","4+")),
    income    = factor(as.character(income),
                       levels = c("$4,000-$5,999","<=$1,999",
                                  "$2,000-$3,999","$6,000-$7,999",
                                  "$8,000-$9,999","$10,000-$14,999",
                                  ">=$15,000")),
    educ      = factor(educ,      levels = c("<HS","Some HS","HS grad","College+")),
    marital   = factor(marital,   levels = c("Married/CL","Single","Wid/Div/Sep")),
    race      = factor(race,      levels = c("White","Black","Puerto Rican","Other")),
    site      = factor(site)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex, mvm2, ai), as.factor)) %>%
  mutate(ih_obs = as.integer(!is.na(ih1) & exclusion == 0))

imp_data <- weighted %>%
  filter(.imp > 0, ih_obs == 1) %>%
  arrange(MOMID) %>%
  split(.$.imp)

# ------------------------------------------------------------------------------
# RUBIN'S RULES + GENERIC RUNNER (as 10)
# ------------------------------------------------------------------------------
pool_gee <- function(fits) {
  m        <- length(fits)
  coef_mat <- do.call(rbind, lapply(fits, coef))
  se_mat   <- do.call(rbind, lapply(fits, function(f) sqrt(diag(vcov(f)))))
  Q_bar <- colMeans(coef_mat)
  T_var <- colMeans(se_mat^2) + (1 + 1/m) * apply(coef_mat, 2, var)
  SE    <- sqrt(T_var)
  tibble(term = names(Q_bar), log_rr = Q_bar, rr = exp(Q_bar),
         lci = exp(Q_bar - 1.96 * SE), uci = exp(Q_bar + 1.96 * SE), n_imps = m)
}

exposures <- c(mvm2 = "MVM any vs none", mvm = "MVM score (continuous)",
               ai = "AI any vs none")
exp_term  <- c(mvm2 = "mvm21", mvm = "mvm", ai = "ai1")

run <- function(exposure, covs, outcome = "ih1", use_ipw = FALSE,
                subset_fun = identity) {
  rhs <- paste(c(exposure, covs), collapse = " + ")
  f   <- as.formula(paste(outcome, "~", rhs))
  fits <- mclapply(imp_data, function(d) {
    d <- subset_fun(d) %>% mutate(across(where(is.factor), droplevels))
    d$.w <- if (use_ipw) d$ipw else rep(1, nrow(d))
    tryCatch(geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                    corstr = "exchangeable", weights = .w),
             error = function(e) e)
  }, mc.cores = n_cores)
  bad <- sapply(fits, function(x) !inherits(x, "geeglm"))
  if (any(bad)) message("  ", sum(bad), " fit(s) failed: ", rhs)
  d1 <- subset_fun(imp_data[[1]])
  cases_1 <- sum(d1[[outcome]] == 1)
  obs_1   <- nrow(d1)
  pool_gee(fits[!bad]) %>%
    filter(term == exp_term[[exposure]]) %>%
    mutate(exposure = exposures[[exposure]], outcome = outcome,
           n_cases = cases_1, n_obs = obs_1,
           rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci))
}

# ------------------------------------------------------------------------------
# COVARIATES
# ------------------------------------------------------------------------------
cov_labels <- c(age = "Maternal age", bmi = "Pre-pregnancy BMI",
                race = "Race and ethnicity", educ = "Education",
                income = "Income", marital = "Marital status",
                smoking = "Smoking", parity = "Parity",
                dm = "Pre-pregnancy diabetes", chronic_htn = "Chronic hypertension",
                infant_sex = "Infant sex", site = "Study site")
demog   <- c("age", "race", "educ", "income", "marital")
health  <- c("bmi", "smoking", "dm", "chronic_htn")
repro   <- c("parity", "infant_sex")
full    <- c(demog, health, repro)

# ------------------------------------------------------------------------------
# A1. CRUDE + EACH COVARIATE INDIVIDUALLY
# ------------------------------------------------------------------------------
cat("A1: crude + each covariate individually\n")
seq_ind <- map_dfr(names(exposures), function(e) {
  cat(" ", exposures[[e]], "\n")
  crude <- run(e, character(0)) %>% mutate(step = "Crude")
  singles <- map_dfr(names(cov_labels), function(v)
    run(e, v) %>% mutate(step = paste("Crude +", cov_labels[[v]])))
  bind_rows(crude, singles) %>%
    mutate(pct_change_vs_crude =
             sprintf("%+.1f%%", 100 * (log_rr - log_rr[1]) / log_rr[1]))
}) %>%
  select(exposure, step, n_cases, n_obs, rr_ci, pct_change_vs_crude, rr, lci, uci)

# ------------------------------------------------------------------------------
# A2. CUMULATIVE BLOCKS + PRIMARY AND SECONDARY
# ------------------------------------------------------------------------------
cat("A2: cumulative adjustment\n")
steps <- list(
  "1. Crude"                                   = character(0),
  "2. + Demographics (age, race, education, income, marital)" = demog,
  "3. + Health (BMI, smoking, diabetes, chronic HTN)"         = c(demog, health),
  "4. + Reproductive (parity, infant sex) = full, no site"    = full,
  "5. + Study site = full with site"                          = c(full, "site")
)
seq_cum <- map_dfr(names(exposures), function(e) {
  cat(" ", exposures[[e]], "\n")
  bind_rows(
    imap_dfr(steps, function(covs, s) run(e, covs) %>% mutate(step = s)),
    run(e, full, use_ipw = TRUE) %>%
      mutate(step = "SECONDARY: full, no site, + IPW"),
    run(e, c(full, "site"), use_ipw = TRUE) %>%
      mutate(step = "PRIMARY: full + site + IPW")
  ) %>%
    mutate(pct_change_vs_crude =
             sprintf("%+.1f%%", 100 * (log_rr - log_rr[1]) / log_rr[1]))
}) %>%
  select(exposure, step, n_cases, n_obs, rr_ci, pct_change_vs_crude, rr, lci, uci)

# ------------------------------------------------------------------------------
# B. BOSTON VS NOT BOSTON (fully adjusted + IPW)
# ------------------------------------------------------------------------------
cat("B: Boston vs not Boston\n")
boston     <- function(d) filter(d, site == "5")
not_boston <- function(d) filter(d, site != "5")
strat <- map_dfr(c("ih1", "ih2"), function(o) map_dfr(names(exposures), function(e) {
  cat(" ", o, exposures[[e]], "\n")
  bind_rows(
    run(e, full, o, TRUE, boston)             %>% mutate(stratum = "Boston (site 5)"),
    run(e, c(full, "site"), o, TRUE, not_boston) %>% mutate(stratum = "Other 11 sites (+ site)"),
    run(e, c(full, "site"), o, TRUE)          %>% mutate(stratum = "All sites (primary)")
  )
})) %>%
  mutate(outcome = recode(outcome, ih1 = "Definite IH", ih2 = "Suspect = present")) %>%
  select(outcome, exposure, stratum, n_cases, n_obs, rr_ci, rr, lci, uci)

# ------------------------------------------------------------------------------
# PRINT + SAVE
# ------------------------------------------------------------------------------
cat("\n==================== A1. INDIVIDUAL COVARIATES ====================\n")
print(seq_ind %>% select(exposure, step, rr_ci, pct_change_vs_crude), n = Inf, width = Inf)
cat("\n==================== A2. CUMULATIVE ====================\n")
print(seq_cum %>% select(exposure, step, n_cases, n_obs, rr_ci, pct_change_vs_crude),
      n = Inf, width = Inf)
cat("\n==================== B. BOSTON VS NOT BOSTON ====================\n")
print(strat %>% select(outcome, exposure, stratum, n_cases, n_obs, rr_ci), n = Inf, width = Inf)

write_csv(seq_ind, file.path(out_dir, "seq_individual.csv"))
write_csv(seq_cum, file.path(out_dir, "seq_cumulative.csv"))
write_csv(strat,   file.path(out_dir, "boston_stratified.csv"))
cat("\nSaved: outputs/seq_individual.csv, seq_cumulative.csv, boston_stratified.csv\n")
