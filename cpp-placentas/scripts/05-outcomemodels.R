# ==============================================================================
# 05-models.R
# CPP Placental Pathology → Infantile Hemangioma
# Updated May 2026
#
# GEE Poisson with robust sandwich SE, clustered on MOMID
# Pooled across 20 imputations using Rubin's rules
# Primary/sensitivity: exchangeable correlation
# Boston sensitivity: independence correlation (faster, same point estimates)
# droplevels() applied after site filtering to avoid unused factor level error
# ==============================================================================

library(tidyverse)
library(geepack)

setwd("/Users/wongjj/Downloads")

weighted <- readRDS("weighted_long.RDS")

# ------------------------------------------------------------------------------
# FACTOR CODING
# ------------------------------------------------------------------------------
weighted <- weighted %>%
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
  mutate(across(c(dm, chronic_htn, infant_sex, chorangioma,
                  mvm2, mvm_villous, mvm_vascular, ai), as.factor)) %>%
  mutate(ih_obs = as.integer(!is.na(ih1) & exclusion == 0))

cat("Weighted long dims:", nrow(weighted), "x", ncol(weighted), "\n")
cat("Imputations:", max(weighted$.imp), "\n")
cat("Observed outcome rows per imp:",
    sum(weighted$ih_obs[weighted$.imp==1]==1), "\n")

# ------------------------------------------------------------------------------
# RUBIN'S RULES POOLING FOR GEE
# ------------------------------------------------------------------------------
pool_gee <- function(fits) {
  m        <- length(fits)
  coef_mat <- do.call(rbind, lapply(fits, coef))
  se_mat   <- do.call(rbind, lapply(fits, function(f) sqrt(diag(vcov(f)))))
  Q_bar <- colMeans(coef_mat)
  U_bar <- colMeans(se_mat^2)
  B     <- apply(coef_mat, 2, var)
  T_var <- U_bar + (1 + 1/m) * B
  SE    <- sqrt(T_var)
  z     <- Q_bar / SE
  pval  <- 2 * pnorm(abs(z), lower.tail = FALSE)
  tibble(
    term    = names(Q_bar),
    log_rr  = Q_bar,
    se      = SE,
    rr      = exp(Q_bar),
    lci     = exp(Q_bar - 1.96 * SE),
    uci     = exp(Q_bar + 1.96 * SE),
    p_value = pval
  )
}

# ------------------------------------------------------------------------------
# COVARIATE SETS
# ------------------------------------------------------------------------------
confounders <- "age + bmi + race + educ + income + marital +
                smoking + parity + dm + chronic_htn + infant_sex + site"

confounders_nosite <- "age + bmi + race + educ + income + marital +
                       smoking + parity + dm + chronic_htn + infant_sex"

cov_vars <- c("age", "bmi", "race", "educ", "income", "marital",
              "smoking", "parity", "dm", "chronic_htn", "infant_sex", "site")

cov_vars_nosite <- c("age", "bmi", "race", "educ", "income", "marital",
                     "smoking", "parity", "dm", "chronic_htn", "infant_sex")

# ------------------------------------------------------------------------------
# GENERIC MODEL RUNNER
# corstr = "exchangeable" for main models
# corstr = "independence" for Boston sensitivity (avoids hanging)
# droplevels() after site filtering prevents unused factor level error
# ------------------------------------------------------------------------------
run_gee_generic <- function(exposure, label, outcome,
                            use_ipw     = FALSE,
                            conf_str    = confounders,
                            covs        = cov_vars,
                            site_filter = NULL,
                            corstr      = "exchangeable") {
  m    <- max(weighted$.imp)
  fits <- vector("list", m)
  exp_vars <- trimws(unlist(strsplit(exposure, "\\*")))
  exp_vars <- exp_vars[exp_vars %in% names(weighted)]
  
  for (i in seq_len(m)) {
    d <- weighted %>% filter(.imp == i, ih_obs == 1)
    
    if (!is.null(site_filter)) {
      if (site_filter == "boston_only") d <- d %>% filter(site == 5)
      if (site_filter == "excl_boston") d <- d %>% filter(site != 5)
    }
    
    d <- d %>%
      filter(if_all(all_of(c(covs, exp_vars)), ~ !is.na(.))) %>%
      mutate(across(where(is.factor), droplevels)) %>%
      arrange(MOMID)
    
    f <- as.formula(paste(outcome, "~", exposure, "+", conf_str))
    w <- if (use_ipw) d$ipw else rep(1, nrow(d))
    
    fits[[i]] <- tryCatch(
      geeglm(f, data = d, family = poisson(link = "log"),
             id = MOMID, corstr = corstr, weights = w),
      error = function(e) {
        message("  Error imp ", i, ": ", e$message); NULL
      }
    )
  }
  
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) { warning("All fits failed: ", label); return(NULL) }
  
  pool_gee(fits) %>%
    filter(
      term != "(Intercept)",
      !grepl("^site|^race|^educ|^income|^marital|^smoking|^parity|^age$|^bmi$|^dm|^chronic_htn|^infant_sex",
             term)
    ) %>%
    mutate(
      exposure = label,
      outcome  = outcome,
      rr_ci    = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci),
      p_fmt    = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
    )
}

# ------------------------------------------------------------------------------
# WRAPPERS
# ------------------------------------------------------------------------------
run_gee <- function(e, l, o, ipw = FALSE)
  run_gee_generic(e, l, o, ipw,
                  conf_str = confounders, covs = cov_vars,
                  site_filter = NULL, corstr = "exchangeable") %>%
  mutate(model = ifelse(ipw, "adjusted+IPW", "adjusted"))

run_gee_nosite <- function(e, l, o, ipw = FALSE)
  run_gee_generic(e, l, o, ipw,
                  conf_str = confounders_nosite, covs = cov_vars_nosite,
                  site_filter = NULL, corstr = "exchangeable") %>%
  mutate(model = ifelse(ipw, "adjusted+IPW (no site)", "adjusted (no site)"))

run_gee_boston <- function(e, l, o, ipw = FALSE)
  run_gee_generic(e, l, o, ipw,
                  conf_str = confounders_nosite, covs = cov_vars_nosite,
                  site_filter = "boston_only", corstr = "independence") %>%
  mutate(model = ifelse(ipw, "Boston only + IPW", "Boston only"))

run_gee_noboston <- function(e, l, o, ipw = FALSE)
  run_gee_generic(e, l, o, ipw,
                  conf_str = confounders, covs = cov_vars,
                  site_filter = "excl_boston", corstr = "independence") %>%
  mutate(model = ifelse(ipw, "Excl Boston + IPW", "Excl Boston"))

# ------------------------------------------------------------------------------
# EXPOSURES
# ------------------------------------------------------------------------------
exposures <- list(
  list(var = "mvm3",        label = "MVM grade (3-level)"),
  list(var = "mvm",         label = "MVM score (continuous)"),
  list(var = "mvm2",        label = "MVM any vs none"),
  list(var = "ai3",         label = "AI stage (3-level)"),
  list(var = "ai",          label = "AI any vs none"),
  list(var = "mvm2 * ai",   label = "MVM x AI joint"),
  list(var = "chorangioma", label = "Chorangioma")
)

# ------------------------------------------------------------------------------
# RUN ALL MODELS
# ------------------------------------------------------------------------------
cat("\nFitting primary models (ih1, with site)...\n")
results_v1 <- map_dfr(exposures, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(run_gee(e$var, e$label, "ih1", FALSE),
            run_gee(e$var, e$label, "ih1", TRUE))
}) %>% mutate(outcome_label = "V1: suspect -> absent (primary)")

cat("\nFitting sensitivity — ih2...\n")
results_v2 <- map_dfr(exposures, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(run_gee(e$var, e$label, "ih2", FALSE),
            run_gee(e$var, e$label, "ih2", TRUE))
}) %>% mutate(outcome_label = "V2: suspect -> present (sensitivity)")

cat("\nFitting sensitivity — no site...\n")
results_nosite <- map_dfr(exposures, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(run_gee_nosite(e$var, e$label, "ih1", FALSE),
            run_gee_nosite(e$var, e$label, "ih1", TRUE))
}) %>% mutate(outcome_label = "V1: no site adjustment")

cat("\nFitting sensitivity — Boston only (independence)...\n")
results_boston <- map_dfr(exposures, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(run_gee_boston(e$var, e$label, "ih1", FALSE),
            run_gee_boston(e$var, e$label, "ih1", TRUE))
}) %>% mutate(outcome_label = "V1: Boston only (site 5)")

cat("\nFitting sensitivity — exclude Boston (independence)...\n")
results_noboston <- map_dfr(exposures, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(run_gee_noboston(e$var, e$label, "ih1", FALSE),
            run_gee_noboston(e$var, e$label, "ih1", TRUE))
}) %>% mutate(outcome_label = "V1: exclude Boston")

# ------------------------------------------------------------------------------
# COMBINE
# ------------------------------------------------------------------------------
results <- bind_rows(results_v1, results_v2, results_nosite,
                     results_boston, results_noboston) %>%
  select(outcome_label, exposure, model, term,
         rr, lci, uci, rr_ci, p_value, p_fmt)

# ------------------------------------------------------------------------------
# PRINT
# ------------------------------------------------------------------------------
cat("\n======================================================================\n")
cat("PRIMARY — ih1, adjusted + IPW (with site)\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V1: suspect -> absent (primary)",
         model == "adjusted+IPW") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("SENSITIVITY — ih2, adjusted + IPW\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V2: suspect -> present (sensitivity)",
         model == "adjusted+IPW") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("SENSITIVITY — no site, adjusted + IPW\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V1: no site adjustment",
         model == "adjusted+IPW (no site)") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("SENSITIVITY — Boston only, adjusted + IPW\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V1: Boston only (site 5)",
         model == "Boston only + IPW") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("SENSITIVITY — exclude Boston, adjusted + IPW\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V1: exclude Boston",
         model == "Excl Boston + IPW") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("THREE-WAY: Full / Boston only / Exclude Boston\n")
cat("======================================================================\n")
bind_rows(
  results %>%
    filter(outcome_label == "V1: suspect -> absent (primary)",
           model == "adjusted+IPW") %>%
    mutate(model = "Full sample"),
  results %>%
    filter(outcome_label == "V1: Boston only (site 5)",
           model == "Boston only + IPW") %>%
    mutate(model = "Boston only"),
  results %>%
    filter(outcome_label == "V1: exclude Boston",
           model == "Excl Boston + IPW") %>%
    mutate(model = "Excl Boston")
) %>%
  select(model, exposure, term, rr_ci, p_fmt) %>%
  arrange(exposure, term, model) %>%
  print(n = Inf)

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
write_csv(results, "outcome_results.csv")
saveRDS(results,   "outcome_results.RDS")

cat("\nSaved: outcome_results.csv, outcome_results.RDS\n")
cat("Total rows:", nrow(results), "\n")