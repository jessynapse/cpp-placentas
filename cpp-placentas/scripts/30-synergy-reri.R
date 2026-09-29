# ==============================================================================
# 30-synergy-reri.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Synergy between MVM and AI, requested by Linda Ernst ("Is there any
# synergy between acute inflammation and MVM?"). Supplemental table.
#
# Uses the four-level joint exposure from 16 (neither = reference, MVM
# only, AI only, both), as recommended for reporting interaction by Knol
# and VanderWeele (Int J Epidemiol 2012), and adds:
#   Additive scale
#     RERI = RR11 - RR10 - RR01 + 1    (0 = no additive interaction)
#     AP   = RERI / RR11               (proportion of risk in the doubly
#                                       exposed due to interaction)
#     S    = (RR11 - 1) / ((RR10 - 1) + (RR01 - 1))   (1 = none)
#   Multiplicative scale
#     RR11 / (RR10 x RR01)             (1 = none)
# 95% CIs by the delta method (Hosmer and Lemeshow 1992) on the pooled
# coefficients and the pooled variance-covariance matrix. The covariance
# matrix is pooled with Rubin's rules: mean within-imputation covariance
# + (1 + 1/m) x between-imputation covariance of the coefficients. S is
# computed on the log scale. S is unstable when RR10 and RR01 are both
# close to 1, so its CI should be read with caution.
#
# Corrected pipeline (07-10): 40,700 infants, 700 definite IH cases,
# 20 imputations, modified Poisson GEE (exchangeable, clustered on MOMID).
# Models: crude, primary (11 covariates + site + IPW), secondary (11
# covariates + IPW, no site).
#
# Output: outputs/synergy_reri.csv. No p-values.
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
# DATA (as 16)
# ------------------------------------------------------------------------------
imp_data <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  filter(.imp > 0, !is.na(ih1), exclusion == 0) %>%
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
    joint     = factor(case_when(mvm2 == 0 & ai == 0 ~ "Neither",
                                 mvm2 == 1 & ai == 0 ~ "MVM only",
                                 mvm2 == 0 & ai == 1 ~ "AI only",
                                 mvm2 == 1 & ai == 1 ~ "Both"),
                       levels = c("Neither", "MVM only", "AI only", "Both"))
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor)) %>%
  select(.imp, MOMID, ih1, ipw, joint, site, age, bmi, race, educ,
         income, marital, smoking, parity, dm, chronic_htn, infant_sex) %>%
  arrange(MOMID) %>%
  split(.$.imp)
invisible(gc())
stopifnot(nrow(imp_data[[1]]) == 40700, sum(imp_data[[1]]$ih1) == 700)

# ------------------------------------------------------------------------------
# MODELS: pooled coefficients AND pooled covariance matrix
# ------------------------------------------------------------------------------
confounders <- "age + bmi + race + educ + income + marital +
                smoking + parity + dm + chronic_htn + infant_sex"
models <- c("Crude" = NA, "Primary: covariates + site + IPW" = "+ site",
            "Secondary: covariates + IPW, no site" = "")
jt <- c("jointMVM only", "jointAI only", "jointBoth")

pool_joint <- function(extra) {
  rhs <- if (is.na(extra)) "joint" else paste("joint +", confounders, extra)
  f <- as.formula(paste("ih1 ~", rhs))
  use_ipw <- !is.na(extra)
  fit_one <- function(d) {
    d$.w <- if (use_ipw) d$ipw else rep(1, nrow(d))
    tryCatch({
      fit <- geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                    corstr = "exchangeable", weights = .w)
      list(b = coef(fit)[jt], V = vcov(fit)[jt, jt])
    }, error = function(err) err)
  }
  fits <- mclapply(imp_data, fit_one, mc.cores = n_cores)
  redo <- which(!sapply(fits, is.list) | sapply(fits, inherits, "error"))
  for (i in redo) fits[[i]] <- fit_one(imp_data[[i]])
  ok <- sapply(fits, function(x) is.list(x) && !inherits(x, "error"))
  if (any(!ok)) message("  ", sum(!ok), " fit(s) failed: ", rhs)
  fits <- fits[ok]; m <- length(fits)
  B <- do.call(rbind, lapply(fits, `[[`, "b"))
  list(b = colMeans(B),
       V = Reduce(`+`, lapply(fits, `[[`, "V")) / m + (1 + 1/m) * cov(B),
       m = m)
}

# ------------------------------------------------------------------------------
# INTERACTION MEASURES (delta method)
# ------------------------------------------------------------------------------
ci <- function(est, g, V, log_scale = FALSE) {
  se <- sqrt(as.numeric(t(g) %*% V %*% g))
  if (log_scale) c(exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))
  else c(est, est - 1.96 * se, est + 1.96 * se)
}

measures <- function(p) {
  b <- p$b; V <- p$V
  r10 <- exp(b[1]); r01 <- exp(b[2]); r11 <- exp(b[3])
  reri <- r11 - r10 - r01 + 1
  out <- list(
    `RR, MVM only vs neither` = ci(b[1], c(1, 0, 0), V, TRUE),
    `RR, AI only vs neither`  = ci(b[2], c(0, 1, 0), V, TRUE),
    `RR, both vs neither`     = ci(b[3], c(0, 0, 1), V, TRUE),
    `RERI (additive, 0 = none)` =
      ci(reri, c(-r10, -r01, r11), V),
    `Attributable proportion, AP (additive, 0 = none)` =
      ci(reri / r11, c(-r10 / r11, -r01 / r11, (r10 + r01 - 1) / r11), V),
    `Synergy index, S (additive, 1 = none)` =
      ci(log((r11 - 1) / (r10 + r01 - 2)),
         c(-r10 / (r10 + r01 - 2), -r01 / (r10 + r01 - 2), r11 / (r11 - 1)),
         V, TRUE),
    `Ratio of RRs (multiplicative, 1 = none)` =
      ci(b[3] - b[1] - b[2], c(-1, -1, 1), V, TRUE))
  imap_dfr(out, ~ tibble(measure = .y,
                         est = sprintf("%.2f (%.2f, %.2f)", .x[1], .x[2], .x[3])))
}

results <- imap_dfr(models, function(extra, lab) {
  cat(" ", lab, "\n")
  p <- pool_joint(extra)
  measures(p) %>% mutate(model = lab, n_imps = p$m)
})

counts <- imp_data[[1]] %>% group_by(joint) %>%
  summarise(cases = sum(ih1 == 1), obs = n(), .groups = "drop") %>%
  mutate(ih = sprintf("%.1f%% (%s/%s)", 100 * cases / obs,
                      format(cases, big.mark = ","), format(obs, big.mark = ",")))

table_out <- results %>%
  select(measure, model, est) %>%
  pivot_wider(names_from = model, values_from = est) %>%
  mutate(n_cases = 700, n_obs = 40700)

cat("\n==================== IH BY JOINT GROUP ====================\n")
print(counts, width = Inf)
cat("\n==================== MVM x AI INTERACTION (definite IH) ====================\n")
print(table_out, width = Inf)

write_csv(table_out, file.path(out_dir, "synergy_reri.csv"))
cat("\nSaved: outputs/synergy_reri.csv\n")
