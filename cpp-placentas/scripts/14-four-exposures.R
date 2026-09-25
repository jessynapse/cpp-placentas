# ==============================================================================
# 14-four-exposures.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Main results table for four exposures, definite IH (ih1):
#   1. MVM any vs none (mvm2)
#   2. MVM score, per 1 point (mvm, 0 to 8)
#   3. AI any vs none (ai)
#   4. AI score, per 1 compartment (ai_count, 0 to 7)  NEW
#
# AI score = number of placental compartments with any neutrophilic
# infiltration (grade > 0) in path2all.sas7bdat:
#   maternal: amnion membrane roll (PA02_34), chorion membrane roll
#             (PA02_35), amnion placental surface (PA02_36), chorion
#             placental surface (PA02_37)
#   fetal:    umbilical vein (PA02_25), umbilical artery (PA02_26),
#             fetal surface vessels (PA02_38)
# Missing compartment grades count as 0 (not involved). AI_DI equals
# "any of these 7 compartments involved" for 100% of placentas, so
# ai_count >= 1 is identical to ai == 1. This is checked below.
# The score parallels the MVM score (a count of lesions present).
#
# Models (corrected pipeline 07-10, 40,700 infants, 700 cases):
#   crude, secondary (11 covariates + IPW), primary (+ site + IPW)
#   Weights from 09 (denominator includes mvm2, ai, covariates, site)
#
# Outputs (outputs/):
#   four_exposures_results.csv   estimates with case and observation counts
#   ai_score_distribution.csv    IH by AI score (and MVM score for comparison)
# No p-values.
# ==============================================================================

library(tidyverse)
library(haven)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"   # folder with path2all.sas7bdat and weighted_long_corrected.RDS
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE)
n_cores <- min(4, detectCores())   # set to 1 on Windows

# ------------------------------------------------------------------------------
# BUILD AI SCORE FROM COMPARTMENT DATA
# ------------------------------------------------------------------------------
ai_comps <- c("PA02_34", "PA02_35", "PA02_36", "PA02_37",
              "PA02_25", "PA02_26", "PA02_38")

path2 <- read_sas(file.path(in_dir, "path2all.sas7bdat"),
                  col_select = all_of(c("MOMID", "PREGID", "CHILDID", ai_comps))) %>%
  zap_labels() %>%
  mutate(ai_count = rowSums(across(all_of(ai_comps), ~ !is.na(.x) & .x > 0))) %>%
  select(MOMID, PREGID, CHILDID, ai_count)

# ------------------------------------------------------------------------------
# DATA (as 10) + AI SCORE
# ------------------------------------------------------------------------------
weighted <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  left_join(path2, by = c("MOMID", "PREGID", "CHILDID")) %>%
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
    site      = factor(site),
    ih_obs    = as.integer(!is.na(ih1) & exclusion == 0),
    mvm2_num  = as.numeric(as.character(mvm2)),
    ai_num    = as.numeric(as.character(ai))
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor))

# Check: AI score reproduces binary AI exactly
chk <- weighted %>% filter(.imp == 0)
stopifnot(!any(is.na(chk$ai_count)))
stopifnot(all((chk$ai_count >= 1) == (chk$ai_num == 1)))
cat("Check passed: ai_count >= 1 identical to AI_DI for all",
    nrow(chk), "pregnancies\n")

imp_data <- weighted %>%
  filter(.imp > 0, ih_obs == 1) %>%
  arrange(MOMID) %>%
  split(.$.imp)

# ------------------------------------------------------------------------------
# RUBIN'S RULES + RUNNER (as 10)
# ------------------------------------------------------------------------------
pool_gee <- function(fits) {
  m        <- length(fits)
  coef_mat <- do.call(rbind, lapply(fits, coef))
  se_mat   <- do.call(rbind, lapply(fits, function(f) sqrt(diag(vcov(f)))))
  Q_bar <- colMeans(coef_mat)
  T_var <- colMeans(se_mat^2) + (1 + 1/m) * apply(coef_mat, 2, var)
  SE    <- sqrt(T_var)
  tibble(term = names(Q_bar), rr = exp(Q_bar),
         lci = exp(Q_bar - 1.96 * SE), uci = exp(Q_bar + 1.96 * SE), n_imps = m)
}

confounders <- "age + bmi + race + educ + income + marital +
                smoking + parity + dm + chronic_htn + infant_sex"

run <- function(exposure, model) {
  rhs <- switch(model,
    "Crude"                      = exposure,
    "Secondary: no site + IPW"   = paste(exposure, "+", confounders),
    "Primary: + site + IPW"      = paste(exposure, "+", confounders, "+ site"))
  f <- as.formula(paste("ih1 ~", rhs))
  use_ipw <- model != "Crude"
  fits <- mclapply(imp_data, function(d) {
    d$.w <- if (use_ipw) d$ipw else rep(1, nrow(d))
    tryCatch(geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                    corstr = "exchangeable", weights = .w),
             error = function(e) e)
  }, mc.cores = n_cores)
  bad <- sapply(fits, function(x) !inherits(x, "geeglm"))
  if (any(bad)) message("  ", sum(bad), " fit(s) failed: ", exposure, " ", model)
  d1 <- imp_data[[1]]
  cases_1 <- sum(d1$ih1 == 1)
  obs_1   <- nrow(d1)
  pool_gee(fits[!bad]) %>%
    filter(term == exposure) %>%
    mutate(model = model, n_cases = cases_1, n_obs = obs_1,
           rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci))
}

exposures <- c(mvm2_num = "MVM any vs none",
               mvm      = "MVM score, per 1 point (0 to 8)",
               ai_num   = "AI any vs none",
               ai_count = "AI score, per 1 compartment (0 to 7)")
models <- c("Crude", "Secondary: no site + IPW", "Primary: + site + IPW")

results <- map_dfr(names(exposures), function(e) map_dfr(models, function(m) {
  cat(" ", exposures[[e]], "/", m, "\n")
  run(e, m) %>% mutate(exposure = exposures[[e]])
})) %>%
  select(exposure, model, n_cases, n_obs, rr_ci, rr, lci, uci, n_imps)

# ------------------------------------------------------------------------------
# DISTRIBUTION: IH BY SCORE (observed outcome, never imputed)
# ------------------------------------------------------------------------------
fmt <- function(x, n) sprintf("%.1f%% (%s/%s)", 100 * x / n,
                              format(x, big.mark = ",", trim = TRUE),
                              format(n, big.mark = ",", trim = TRUE))
dist <- bind_rows(
  imp_data[[1]] %>% group_by(score = ai_count) %>%
    summarise(n_cases = sum(ih1 == 1), n_obs = n(), .groups = "drop") %>%
    mutate(exposure = "AI score (compartments)"),
  imp_data[[1]] %>% group_by(score = mvm) %>%
    summarise(n_cases = sum(ih1 == 1), n_obs = n(), .groups = "drop") %>%
    mutate(exposure = "MVM score")
) %>%
  mutate(ih_pct = fmt(n_cases, n_obs)) %>%
  select(exposure, score, n_cases, n_obs, ih_pct)

# ------------------------------------------------------------------------------
# PRINT + SAVE
# ------------------------------------------------------------------------------
cat("\n==================== IH BY SCORE ====================\n")
print(dist, n = Inf, width = Inf)
cat("\n==================== FOUR EXPOSURES, DEFINITE IH ====================\n")
print(results %>% select(exposure, model, n_cases, n_obs, rr_ci), n = Inf, width = Inf)

write_csv(results, file.path(out_dir, "four_exposures_results.csv"))
write_csv(dist,    file.path(out_dir, "ai_score_distribution.csv"))
cat("\nSaved: outputs/four_exposures_results.csv, ai_score_distribution.csv\n")
