# ==============================================================================
# 04-ipw.R
# CPP Placental Pathology → Infantile Hemangioma
# Updated May 2026
#
# Single IPW model for non-survival before 1-year examination
# LTFU confirmed non-differential by exposure (Ellen Francis) — not weighted
# Stabilized weights: numerator = intercept only, denominator = full covariates
# Applied within each of 20 imputed datasets
# ==============================================================================

library(tidyverse)
library(mice)

setwd("/Users/wongjj/Downloads")

long <- read_csv("imputed_long.csv", show_col_types = FALSE)

# ------------------------------------------------------------------------------
# FACTOR CODING
# ------------------------------------------------------------------------------
long <- long %>%
  mutate(
    mvm3      = factor(mvm3,      levels = c("None","Low","High")),
    ai3       = factor(ai3,       levels = c("None","Low","High")),
    smoking   = factor(smoking,   levels = c("Non-smoker","1-19/day","20+/day")),
    parity    = factor(parity,    levels = c("0","1","2","3","4+")),
    plurality = factor(plurality, levels = c("Singleton","Multiple")),
    income    = factor(income,    levels = c("<=$1,999","$2,000-$3,999",
                                             "$4,000-$5,999","$6,000-$7,999",
                                             "$8,000-$9,999","$10,000-$14,999",
                                             ">=$15,000")),
    educ      = factor(educ,      levels = c("<HS","Some HS","HS grad","College+")),
    marital   = factor(marital,   levels = c("Married/CL","Single","Wid/Div/Sep")),
    race      = factor(race,      levels = c("White","Black","Puerto Rican","Other")),
    site      = factor(site)
  ) %>%
  mutate(across(c(dm, infant_sex, chorangioma, mvm2,
                  mvm_villous, mvm_vascular, ai), as.factor))

# ------------------------------------------------------------------------------
# OUTCOME OBSERVATION INDICATOR
# 1 = survived to 1y AND has observed IH outcome
# 0 = non-survivor (exclusion == 1) OR LTFU (ih1 NA among survivors)
# exclusion is numeric 0/1 — not string
# ------------------------------------------------------------------------------
long <- long %>%
  mutate(ih_obs = as.integer(!is.na(ih1) & exclusion == 0))

cat("ih_obs distribution (should be same across all imputations):\n")
long %>%
  filter(.imp %in% c(0, 1, 20)) %>%
  group_by(.imp) %>%
  summarise(
    n_observed   = sum(ih_obs == 1),
    n_missing    = sum(ih_obs == 0),
    pct_observed = round(100 * mean(ih_obs == 1), 1),
    .groups = "drop"
  ) %>%
  print()

# ------------------------------------------------------------------------------
# IPW MODEL
# Denominator: P(observed | exposure + covariates)
# Numerator:   P(observed | 1) — stabilized weights
# Weight = numerator / denominator
# newdata = d forces prediction on all rows — returns NA for missing covariates
# ------------------------------------------------------------------------------
ipw_covs <- "age + bmi + race + educ + income + marital +
             smoking + parity + dm + infant_sex + site"

calc_weights <- function(data, imp_num) {
  d <- data %>% filter(.imp == imp_num)
  
  f_denom <- as.formula(paste("ih_obs ~ mvm3 + ai3 +", ipw_covs))
  f_numer <- as.formula("ih_obs ~ 1")
  
  m_denom <- glm(f_denom, data = d, family = binomial(link = "logit"),
                 na.action = na.exclude)
  m_numer <- glm(f_numer, data = d, family = binomial(link = "logit"),
                 na.action = na.exclude)
  
  d$p_obs_denom <- predict(m_denom, newdata = d, type = "response")
  d$p_obs_numer <- predict(m_numer, newdata = d, type = "response")
  d$ipw         <- d$p_obs_numer / d$p_obs_denom
  
  d
}

cat("\nComputing IPW for", max(long$.imp), "imputations...\n")
out <- vector("list", max(long$.imp))
for (i in seq_len(max(long$.imp))) {
  cat("  imputation", i, "\n")
  out[[i]] <- calc_weights(long, i)
}

# Add NA weights for imp 0 (original — not used in outcome models)
orig <- long %>%
  filter(.imp == 0) %>%
  mutate(p_obs_denom = NA_real_,
         p_obs_numer = NA_real_,
         ipw         = NA_real_)

weighted <- bind_rows(orig, bind_rows(out))

# ------------------------------------------------------------------------------
# WEIGHT DIAGNOSTICS
# ------------------------------------------------------------------------------
diag <- weighted %>% filter(.imp > 0)

cat("\n======================================================================\n")
cat("WEIGHT DIAGNOSTICS\n")
cat("======================================================================\n")

cat("\nMean IPW (stabilized — should be ~1.0):",
    round(mean(diag$ipw, na.rm = TRUE), 3), "\n")

summary_weights <- function(x, label) {
  q <- quantile(x, probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99),
                na.rm = TRUE)
  cat(sprintf("\n  %s\n", label))
  cat(sprintf("    Mean: %.3f  Median: %.3f  Max: %.3f\n",
              mean(x, na.rm=TRUE), q["50%"], max(x, na.rm=TRUE)))
  cat(sprintf("    p1: %.3f  p5: %.3f  p95: %.3f  p99: %.3f\n",
              q["1%"], q["5%"], q["95%"], q["99%"]))
  cat(sprintf("    N > 10: %d  N > 20: %d\n",
              sum(x > 10, na.rm=TRUE), sum(x > 20, na.rm=TRUE)))
}

summary_weights(diag$ipw, "All subjects")
summary_weights(diag$ipw[diag$ih_obs == 1], "Observed IH outcome only")

cat("\nMean IPW by imputation (should be consistent ~1.0):\n")
diag %>%
  group_by(.imp) %>%
  summarise(
    mean_ipw = round(mean(ipw, na.rm=TRUE), 3),
    max_ipw  = round(max(ipw, na.rm=TRUE), 2),
    n_gt10   = sum(ipw > 10, na.rm=TRUE),
    .groups  = "drop"
  ) %>%
  print(n = 20)

# ------------------------------------------------------------------------------
# DIAGNOSTIC PLOTS
# ------------------------------------------------------------------------------
pdf("ipw_diagnostics.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))
hist(diag$ipw,
     main = "IPW — all subjects", xlab = "IPW", breaks = 50)
hist(diag$ipw[diag$ih_obs == 1],
     main = "IPW — observed IH outcome only", xlab = "IPW", breaks = 50)
hist(log(diag$ipw + 0.01),
     main = "log(IPW) — all subjects", xlab = "log(IPW)", breaks = 50)
hist(log(diag$ipw[diag$ih_obs == 1] + 0.01),
     main = "log(IPW) — observed IH only", xlab = "log(IPW)", breaks = 50)
dev.off()

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
saveRDS(weighted, "weighted_long.RDS")
write_csv(weighted, "weighted_long.csv", na = "")

cat("\nSaved: weighted_long.RDS, weighted_long.csv\n")
cat("       ipw_diagnostics.pdf\n")
cat("Rows:", nrow(weighted),
    "| Imputations:", max(weighted$.imp), "\n")