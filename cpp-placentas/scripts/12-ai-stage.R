# ==============================================================================
# 12-ai-stage.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Dose-response for acute inflammation (AI) stage, replacing binary AI_DI
# with AI_3cat (0 none, 1 low stage, 2 high stage). Requested by Alexa
# Freedman (comment 6: "AI stage ... might expect to see a dose response").
#
# Same methods as the corrected primary analysis (07-10):
#   - Corrected sample (45,268 singleton pregnancies, complete case on both
#     exposures), 20 imputations from 08
#   - Modified Poisson GEE, exchangeable correlation, clustered on MOMID,
#     pooled with Rubin's rules
#   - Stabilised IPW for an observed IH outcome, REBUILT with ai3 in the
#     denominator model in place of binary ai (mvm2 and all covariates
#     and site as in 09). Numerator intercept only.
#
# Models, for both ih1 (definite) and ih2 (suspect = present):
#   crude, adjusted, adjusted + IPW, adjusted + site + IPW
#   Categorical: low vs none, high vs none (none = reference)
#   Trend: ai3 entered as 0/1/2, RR per one-stage increase
#
# CAUTION: the exact rule used to build AI_3cat is not documented. From
# the compartment data it appears to be a Redline-type stage (high = any
# amnion or umbilical artery involvement), which reproduces ~95% of
# placentas. Confirm the definition with Alexa before describing it.
#
# Outputs (outputs/):
#   ai_stage_results.csv     estimates with case and observation counts
#   ai_stage_counts.csv      IH cases and observations per stage
#   ai_stage_weights.csv     weight diagnostics
# No p-values.
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"   # folder with imputed_long_corrected.csv
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE)
n_cores <- min(4, detectCores())   # set to 1 on Windows

# ------------------------------------------------------------------------------
# DATA + FACTOR CODING (as 09 and 10)
# ------------------------------------------------------------------------------
long <- read_csv(file.path(in_dir, "imputed_long_corrected.csv"),
                 show_col_types = FALSE) %>%
  mutate(
    ai3       = factor(ai3,       levels = c("None","Low","High")),
    smoking   = factor(smoking,   levels = c("Non-smoker","1-19/day","20+/day")),
    parity    = factor(parity,    levels = c("0","1","2","3","4+")),
    income    = factor(income,    levels = c("$4,000-$5,999","<=$1,999",
                                             "$2,000-$3,999","$6,000-$7,999",
                                             "$8,000-$9,999","$10,000-$14,999",
                                             ">=$15,000")),
    educ      = factor(educ,      levels = c("<HS","Some HS","HS grad","College+")),
    marital   = factor(marital,   levels = c("Married/CL","Single","Wid/Div/Sep")),
    race      = factor(race,      levels = c("White","Black","Puerto Rican","Other")),
    site      = factor(site)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex, mvm2), as.factor)) %>%
  mutate(ih_obs     = as.integer(!is.na(ih1) & exclusion == 0),
         ai3_trend  = as.numeric(ai3) - 1)   # 0 none, 1 low, 2 high

stopifnot(!any(is.na(long$ai3)))

# ------------------------------------------------------------------------------
# IPW REBUILT WITH ai3 IN THE DENOMINATOR
# ------------------------------------------------------------------------------
ipw_covs <- "age + bmi + race + educ + income + marital +
             smoking + parity + dm + chronic_htn + infant_sex + site"
f_denom  <- as.formula(paste("ih_obs ~ mvm2 + ai3 +", ipw_covs))

imps <- long %>% filter(.imp > 0) %>% split(.$.imp)
imps <- lapply(imps, function(d) {
  p_den <- predict(glm(f_denom, data = d, family = binomial), type = "response")
  p_num <- mean(d$ih_obs)
  d$ipw_ai3 <- p_num / p_den
  d
})

w_all <- unlist(lapply(imps, `[[`, "ipw_ai3"))
w_obs <- unlist(lapply(imps, function(d) d$ipw_ai3[d$ih_obs == 1]))
weights_diag <- tibble(
  subset = c("All subjects", "Observed IH outcome"),
  mean   = c(mean(w_all), mean(w_obs)),
  min    = c(min(w_all),  min(w_obs)),
  p1     = c(quantile(w_all, .01), quantile(w_obs, .01)),
  p99    = c(quantile(w_all, .99), quantile(w_obs, .99)),
  max    = c(max(w_all),  max(w_obs))
) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
cat("Weight diagnostics (ai3 in denominator):\n"); print(weights_diag)

imp_data <- lapply(imps, function(d) d %>% filter(ih_obs == 1) %>% arrange(MOMID))

# ------------------------------------------------------------------------------
# RUBIN'S RULES (as 10)
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

run <- function(exposure, outcome, model) {
  rhs <- switch(model,
    "crude"             = exposure,
    "adjusted"          = paste(exposure, "+", confounders),
    "adjusted+IPW"      = paste(exposure, "+", confounders),
    "adjusted+site+IPW" = paste(exposure, "+", confounders, "+ site"))
  f <- as.formula(paste(outcome, "~", rhs))
  use_ipw <- grepl("IPW", model)

  fits <- mclapply(imp_data, function(d) {
    d$.w <- if (use_ipw) d$ipw_ai3 else rep(1, nrow(d))
    tryCatch(geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                    corstr = "exchangeable", weights = .w),
             error = function(e) e)
  }, mc.cores = n_cores)
  bad <- sapply(fits, function(x) !inherits(x, "geeglm"))
  if (any(bad)) message("  ", sum(bad), " fit(s) failed: ", exposure, " ", outcome, " ", model)

  d1 <- imp_data[[1]]
  cases_1 <- sum(d1[[outcome]] == 1)
  obs_1   <- nrow(d1)
  pool_gee(fits[!bad]) %>%
    filter(grepl("^ai3", term)) %>%
    mutate(outcome = outcome, model = model,
           n_cases = cases_1, n_obs = obs_1,
           rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci))
}

# ------------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------------
models <- c("crude", "adjusted", "adjusted+IPW", "adjusted+site+IPW")
results <- list()
for (o in c("ih1", "ih2")) for (m in models) for (e in c("ai3", "ai3_trend")) {
  cat(" ", o, m, e, "\n")
  results[[length(results) + 1]] <- run(e, o, m)
}
results <- bind_rows(results) %>%
  mutate(
    outcome  = recode(outcome, ih1 = "Definite IH", ih2 = "Suspect = present"),
    contrast = recode(term, ai3Low = "Low stage vs none",
                      ai3High = "High stage vs none",
                      ai3_trend = "Per one-stage increase (trend)"),
    model    = factor(model, levels = models)) %>%
  arrange(outcome, contrast != "Low stage vs none", contrast, model) %>%
  select(outcome, contrast, model, n_cases, n_obs, rr_ci, rr, lci, uci, n_imps)

# Cases and observations per stage (outcome and exposure never imputed)
counts <- imp_data[[1]] %>%
  pivot_longer(c(ih1, ih2), names_to = "outcome", values_to = "ih") %>%
  group_by(outcome = recode(outcome, ih1 = "Definite IH", ih2 = "Suspect = present"),
           ai_stage = ai3) %>%
  summarise(n_cases = sum(ih == 1), n_obs = n(), .groups = "drop") %>%
  mutate(ih_pct = sprintf("%.1f%% (%s/%s)", 100 * n_cases / n_obs,
                          format(n_cases, big.mark = ",", trim = TRUE),
                          format(n_obs, big.mark = ",", trim = TRUE)))

# ------------------------------------------------------------------------------
# PRINT + SAVE
# ------------------------------------------------------------------------------
cat("\n==================== IH BY AI STAGE ====================\n")
print(counts, n = Inf, width = Inf)
cat("\n==================== ESTIMATES ====================\n")
print(results %>% select(outcome, contrast, model, n_cases, n_obs, rr_ci),
      n = Inf, width = Inf)

write_csv(results,      file.path(out_dir, "ai_stage_results.csv"))
write_csv(counts,       file.path(out_dir, "ai_stage_counts.csv"))
write_csv(weights_diag, file.path(out_dir, "ai_stage_weights.csv"))
cat("\nSaved: outputs/ai_stage_results.csv, ai_stage_counts.csv, ai_stage_weights.csv\n")
