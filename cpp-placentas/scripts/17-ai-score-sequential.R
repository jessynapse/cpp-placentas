# ==============================================================================
# 17-ai-score-sequential.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Adds the AI score (number of compartments with neutrophils, 0 to 7, as
# built in 14) to the sequential adjustment and Boston vs not Boston
# analyses from 13, then combines everything into four-exposure tables:
#   MVM any, MVM score (per 1 point), AI any, AI score (per 1 compartment)
#
# Same methods as 13: definite IH, 40,700 infants, 700 cases, 20
# imputations, Rubin's rules, modified Poisson GEE (exchangeable, MOMID).
# Sequential steps are unweighted. Final rows and Boston models use IPW.
#
# Inputs: outputs/seq_individual.csv, seq_cumulative.csv,
#         boston_stratified.csv (from 13)
# Outputs (outputs/):
#   seq_individual_4exp.csv, seq_cumulative_4exp.csv, boston_4exp.csv
# No p-values.
# ==============================================================================

library(tidyverse)
library(haven)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"   # folder with path2all.sas7bdat and weighted_long_corrected.RDS
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
n_cores <- min(4, detectCores())   # set to 1 on Windows

# ------------------------------------------------------------------------------
# AI SCORE (as 14)
# ------------------------------------------------------------------------------
ai_comps <- c("PA02_34", "PA02_35", "PA02_36", "PA02_37",
              "PA02_25", "PA02_26", "PA02_38")
path2 <- read_sas(file.path(in_dir, "path2all.sas7bdat"),
                  col_select = all_of(c("MOMID", "PREGID", "CHILDID", ai_comps))) %>%
  zap_labels() %>%
  mutate(ai_count = rowSums(across(all_of(ai_comps), ~ !is.na(.x) & .x > 0))) %>%
  select(MOMID, PREGID, CHILDID, ai_count)

# ------------------------------------------------------------------------------
# DATA (as 13), trimmed to model columns
# ------------------------------------------------------------------------------
imp_data <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  filter(.imp > 0, !is.na(ih1), exclusion == 0) %>%
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
    site      = factor(site)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor)) %>%
  select(.imp, MOMID, ih1, ipw, ai_count, site, age, bmi, race, educ, income,
         marital, smoking, parity, dm, chronic_htn, infant_sex) %>%
  arrange(MOMID) %>%
  split(.$.imp)
stopifnot(!any(is.na(imp_data[[1]]$ai_count)))
invisible(gc())

# ------------------------------------------------------------------------------
# RUNNER (as 13)
# ------------------------------------------------------------------------------
pool_gee <- function(fits) {
  m        <- length(fits)
  coef_mat <- do.call(rbind, lapply(fits, coef))
  se_mat   <- do.call(rbind, lapply(fits, function(f) sqrt(diag(vcov(f)))))
  Q_bar <- colMeans(coef_mat)
  SE    <- sqrt(colMeans(se_mat^2) + (1 + 1/m) * apply(coef_mat, 2, var))
  tibble(term = names(Q_bar), log_rr = Q_bar, rr = exp(Q_bar),
         lci = exp(Q_bar - 1.96 * SE), uci = exp(Q_bar + 1.96 * SE))
}

run <- function(covs, use_ipw = FALSE, subset_fun = identity) {
  f <- as.formula(paste("ih1 ~", paste(c("ai_count", covs), collapse = " + ")))
  fit_one <- function(d) {
    d <- subset_fun(d) %>% mutate(across(where(is.factor), droplevels))
    d$.w <- if (use_ipw) d$ipw else rep(1, nrow(d))
    tryCatch(geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                    corstr = "exchangeable", weights = .w),
             error = function(err) err)
  }
  fits <- mclapply(imp_data, fit_one, mc.cores = n_cores)
  redo <- which(!sapply(fits, function(x) inherits(x, "geeglm")))
  for (i in redo) fits[[i]] <- fit_one(imp_data[[i]])
  bad <- sapply(fits, function(x) !inherits(x, "geeglm"))
  if (any(bad)) message("  ", sum(bad), " fit(s) failed")
  d1 <- subset_fun(imp_data[[1]])
  cases_1 <- sum(d1$ih1 == 1); obs_1 <- nrow(d1)
  pool_gee(fits[!bad]) %>% filter(term == "ai_count") %>%
    mutate(exposure = "AI score (continuous)", n_cases = cases_1, n_obs = obs_1,
           rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci))
}

cov_labels <- c(age = "Maternal age", bmi = "Pre-pregnancy BMI",
                race = "Race and ethnicity", educ = "Education",
                income = "Income", marital = "Marital status",
                smoking = "Smoking", parity = "Parity",
                dm = "Pre-pregnancy diabetes", chronic_htn = "Chronic hypertension",
                infant_sex = "Infant sex", site = "Study site")
demog  <- c("age", "race", "educ", "income", "marital")
health <- c("bmi", "smoking", "dm", "chronic_htn")
repro  <- c("parity", "infant_sex")
full   <- c(demog, health, repro)

# A1. individual
cat("A1\n")
ind <- bind_rows(run(character(0)) %>% mutate(step = "Crude"),
                 imap_dfr(cov_labels, function(lab, v) run(v) %>%
                            mutate(step = paste("Crude +", lab))))

# A2. cumulative (step labels exactly as in 13)
cat("A2\n")
steps <- list(
  "1. Crude"                                   = character(0),
  "2. + Demographics (age, race, education, income, marital)" = demog,
  "3. + Health (BMI, smoking, diabetes, chronic HTN)"         = c(demog, health),
  "4. + Reproductive (parity, infant sex) = full, no site"    = full,
  "5. + Study site = full with site"                          = c(full, "site"))
cum <- bind_rows(
  imap_dfr(steps, function(covs, s) run(covs) %>% mutate(step = s)),
  run(full, TRUE)          %>% mutate(step = "SECONDARY: full, no site, + IPW"),
  run(c(full, "site"), TRUE) %>% mutate(step = "PRIMARY: full + site + IPW"))

# B. Boston
cat("B\n")
bos <- bind_rows(
  run(full, TRUE, function(d) filter(d, site == "5")) %>% mutate(stratum = "Boston (site 5)"),
  run(c(full, "site"), TRUE, function(d) filter(d, site != "5")) %>%
    mutate(stratum = "Other 11 sites (+ site)"),
  run(c(full, "site"), TRUE) %>% mutate(stratum = "All sites (primary)"))

# ------------------------------------------------------------------------------
# COMBINE WITH 13 (definite IH only) AND SAVE
# ------------------------------------------------------------------------------
exp_order <- c("MVM any vs none", "MVM score (continuous)", "AI any vs none",
               "AI score (continuous)")
wide <- function(d, key) d %>%
  select(all_of(key), exposure, rr_ci) %>%
  mutate(exposure = factor(exposure, exp_order)) %>%
  arrange(exposure) %>%
  pivot_wider(names_from = exposure, values_from = rr_ci)

s13_ind <- read_csv(file.path(out_dir, "seq_individual.csv"), show_col_types = FALSE)
s13_cum <- read_csv(file.path(out_dir, "seq_cumulative.csv"), show_col_types = FALSE)
b13     <- read_csv(file.path(out_dir, "boston_stratified.csv"), show_col_types = FALSE) %>%
  filter(outcome == "Definite IH")

ind4 <- bind_rows(s13_ind %>% select(exposure, step, rr_ci), ind %>% select(exposure, step, rr_ci)) %>%
  mutate(step = factor(step, unique(s13_ind$step))) %>% wide("step") %>% arrange(step)
cum4 <- bind_rows(s13_cum %>% select(exposure, step, rr_ci), cum %>% select(exposure, step, rr_ci)) %>%
  mutate(step = factor(step, unique(s13_cum$step))) %>% wide("step") %>% arrange(step)
bos4 <- bind_rows(b13 %>% select(exposure, stratum, n_cases, n_obs, rr_ci),
                  bos %>% select(exposure, stratum, n_cases, n_obs, rr_ci)) %>%
  mutate(stratum = paste0(stratum, " [", n_cases, " cases / ",
                          format(n_obs, big.mark = ",", trim = TRUE), "]")) %>%
  select(-n_cases, -n_obs) %>% wide("stratum")

cat("\n== Individual ==\n"); print(ind4, n = Inf, width = Inf)
cat("\n== Cumulative ==\n"); print(cum4, n = Inf, width = Inf)
cat("\n== Boston ==\n");     print(bos4, n = Inf, width = Inf)

write_csv(ind4, file.path(out_dir, "seq_individual_4exp.csv"))
write_csv(cum4, file.path(out_dir, "seq_cumulative_4exp.csv"))
write_csv(bos4, file.path(out_dir, "boston_4exp.csv"))
cat("\nSaved: outputs/seq_individual_4exp.csv, seq_cumulative_4exp.csv, boston_4exp.csv\n")
