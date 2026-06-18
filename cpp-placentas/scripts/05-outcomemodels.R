# ==============================================================================
# 05-models.R
# CPP Placental Pathology → Infantile Hemangioma
# Updated June 2026
#
# Per meeting with Stef, Ellen, and Linda (June 2026):
#   - PRIMARY exposure: binary (any vs none) for MVM and AI
#   - PRIMARY model: no site adjustment (site captures pathologist variability)
#   - SENSITIVITY: with site adjustment
#   - SUPPLEMENTAL: 3-level grade, joint MVM x AI, chorangioma, continuous MVM
#   - Chorangioma labeled as "tumor (possible chorangioma)"
#   - Boston sensitivity removed — computationally prohibitive with independence
#     correlation across 20 imputations; Boston finding (AI aRR 1.28, 95% CI
#     1.00, 1.63) from prior run reported narratively in discussion
#   - Crude models added for primary exposures (V1 and V2) for Table 3 and S2
#
# GEE Poisson with robust sandwich SE, clustered on MOMID
# Pooled across 20 imputations using Rubin's rules (manual — GEE not mira-compatible)
# IPW applied to observed outcome subset only (ih_obs == 1)
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
  mutate(
    ih_obs = as.integer(!is.na(ih1) & exclusion == 0),
    mvm_num = as.numeric(mvm),
    mvm2_num = as.numeric(as.character(mvm2)),
    ai_num   = as.numeric(as.character(ai))
  )

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
# PRIMARY = no site (site captures pathologist variability, not true confounding)
# SENSITIVITY = with site
# ------------------------------------------------------------------------------
confounders_primary <- "age + bmi + race + educ + income + marital +
                        smoking + parity + dm + chronic_htn + infant_sex"

confounders_site    <- "age + bmi + race + educ + income + marital +
                        smoking + parity + dm + chronic_htn + infant_sex + site"

cov_vars_primary <- c("age", "bmi", "race", "educ", "income", "marital",
                      "smoking", "parity", "dm", "chronic_htn", "infant_sex")

cov_vars_site    <- c("age", "bmi", "race", "educ", "income", "marital",
                      "smoking", "parity", "dm", "chronic_htn",
                      "infant_sex", "site")

# ------------------------------------------------------------------------------
# GENERIC MODEL RUNNER
# ------------------------------------------------------------------------------
run_gee_generic <- function(exposure, label, outcome,
                            use_ipw     = FALSE,
                            conf_str    = confounders_primary,
                            covs        = cov_vars_primary,
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
# CRUDE MODEL RUNNER
# Uses numeric exposure variables to avoid factor level issues
# ------------------------------------------------------------------------------
run_gee_crude <- function(exposure_num, label, outcome) {
  m    <- max(weighted$.imp)
  fits <- vector("list", m)
  
  for (i in seq_len(m)) {
    d <- weighted %>%
      filter(.imp == i, ih_obs == 1) %>%
      filter(!is.na(.data[[exposure_num]])) %>%
      arrange(MOMID)
    
    f <- as.formula(paste(outcome, "~", exposure_num))
    
    fits[[i]] <- tryCatch(
      geeglm(f, data = d, family = poisson(link = "log"),
             id = MOMID, corstr = "exchangeable"),
      error = function(e) {
        message("  Error imp ", i, ": ", e$message); NULL
      }
    )
  }
  
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) { warning("All crude fits failed: ", label); return(NULL) }
  
  pool_gee(fits) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      exposure = label,
      outcome  = outcome,
      model    = "crude",
      rr_ci    = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci),
      p_fmt    = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
    )
}

# ------------------------------------------------------------------------------
# WRAPPERS
# ------------------------------------------------------------------------------

# Primary — no site
run_gee_primary <- function(e, l, o, ipw = FALSE)
  run_gee_generic(e, l, o, ipw,
                  conf_str = confounders_primary,
                  covs = cov_vars_primary,
                  site_filter = NULL,
                  corstr = "exchangeable") %>%
  mutate(model = ifelse(ipw, "primary+IPW", "primary"))

# Sensitivity — with site
run_gee_site <- function(e, l, o, ipw = FALSE)
  run_gee_generic(e, l, o, ipw,
                  conf_str = confounders_site,
                  covs = cov_vars_site,
                  site_filter = NULL,
                  corstr = "exchangeable") %>%
  mutate(model = ifelse(ipw, "site+IPW", "site adjusted"))

# ------------------------------------------------------------------------------
# EXPOSURE LISTS
# ------------------------------------------------------------------------------

# Primary exposures — binary + continuous
exposures_primary <- list(
  list(var = "mvm2", label = "MVM any vs none"),
  list(var = "mvm",  label = "MVM score (continuous)"),
  list(var = "ai",   label = "AI any vs none")
)

# Crude exposure numeric names — avoid factor issues
exposures_crude <- list(
  list(var = "mvm2_num", label = "MVM any vs none"),
  list(var = "mvm_num",  label = "MVM score (continuous)"),
  list(var = "ai_num",   label = "AI any vs none")
)

# Supplemental exposures — 3-level, joint, chorangioma
exposures_supp <- list(
  list(var = "mvm3",        label = "MVM grade (3-level)"),
  list(var = "ai3",         label = "AI stage (3-level)"),
  list(var = "mvm2 * ai",   label = "MVM x AI joint"),
  list(var = "chorangioma", label = "Tumor (possible chorangioma)")
)

# ------------------------------------------------------------------------------
# RUN CRUDE MODELS — V1 and V2
# ------------------------------------------------------------------------------
cat("\nFitting CRUDE models (ih1)...\n")
results_crude_v1 <- map_dfr(exposures_crude, function(e) {
  cat(" ", e$label, "\n")
  run_gee_crude(e$var, e$label, "ih1")
}) %>% mutate(outcome_label = "V1: crude")

cat("\nFitting CRUDE models (ih2)...\n")
results_crude_v2 <- map_dfr(exposures_crude, function(e) {
  cat(" ", e$label, "\n")
  run_gee_crude(e$var, e$label, "ih2")
}) %>% mutate(outcome_label = "V2: crude")

# ------------------------------------------------------------------------------
# RUN PRIMARY MODELS — binary MVM and AI, no site, ih1
# ------------------------------------------------------------------------------
cat("\nFitting PRIMARY models (ih1, no site)...\n")
results_primary <- map_dfr(exposures_primary, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(
    run_gee_primary(e$var, e$label, "ih1", ipw = FALSE),
    run_gee_primary(e$var, e$label, "ih1", ipw = TRUE)
  )
}) %>% mutate(outcome_label = "V1: primary (no site)")

# ------------------------------------------------------------------------------
# SENSITIVITY 1 — with site adjustment
# ------------------------------------------------------------------------------
cat("\nFitting SENSITIVITY: site adjustment (ih1)...\n")
results_site <- map_dfr(exposures_primary, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(
    run_gee_site(e$var, e$label, "ih1", ipw = FALSE),
    run_gee_site(e$var, e$label, "ih1", ipw = TRUE)
  )
}) %>% mutate(outcome_label = "V1: sensitivity (with site)")

# ------------------------------------------------------------------------------
# SENSITIVITY 2 — V2 outcome (suspect -> present), no site
# ------------------------------------------------------------------------------
cat("\nFitting SENSITIVITY: V2 outcome (suspect -> present)...\n")
results_v2 <- map_dfr(exposures_primary, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(
    run_gee_primary(e$var, e$label, "ih2", ipw = FALSE),
    run_gee_primary(e$var, e$label, "ih2", ipw = TRUE)
  )
}) %>% mutate(outcome_label = "V2: suspect -> present")

# ------------------------------------------------------------------------------
# SUPPLEMENTAL MODELS — 3-level, joint, chorangioma
# ------------------------------------------------------------------------------
cat("\nFitting SUPPLEMENTAL models (ih1, no site)...\n")
results_supp <- map_dfr(exposures_supp, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(
    run_gee_primary(e$var, e$label, "ih1", ipw = FALSE),
    run_gee_primary(e$var, e$label, "ih1", ipw = TRUE)
  )
}) %>% mutate(outcome_label = "V1: supplemental")

cat("\nFitting SUPPLEMENTAL models V2 (suspect -> present)...\n")
results_supp_v2 <- map_dfr(exposures_supp, function(e) {
  cat(" ", e$label, "\n")
  bind_rows(
    run_gee_primary(e$var, e$label, "ih2", ipw = FALSE),
    run_gee_primary(e$var, e$label, "ih2", ipw = TRUE)
  )
}) %>% mutate(outcome_label = "V2: supplemental")

# ------------------------------------------------------------------------------
# COMBINE ALL
# ------------------------------------------------------------------------------
results <- bind_rows(
  results_crude_v1,
  results_crude_v2,
  results_primary,
  results_site,
  results_v2,
  results_supp,
  results_supp_v2
) %>%
  select(outcome_label, exposure, model, term,
         rr, lci, uci, rr_ci, p_value, p_fmt)

# ------------------------------------------------------------------------------
# PRINT KEY RESULTS
# ------------------------------------------------------------------------------
cat("\n======================================================================\n")
cat("CRUDE — ih1\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V1: crude") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("CRUDE — ih2\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V2: crude") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("PRIMARY — ih1, no site, adjusted + IPW\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V1: primary (no site)",
         model == "primary+IPW") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("SENSITIVITY — ih1, with site, adjusted + IPW\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V1: sensitivity (with site)",
         model == "site+IPW") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("SENSITIVITY — V2 (suspect -> present), no site, adjusted + IPW\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V2: suspect -> present",
         model == "primary+IPW") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("SUPPLEMENTAL — 3-level, joint, chorangioma (no site, adjusted + IPW)\n")
cat("======================================================================\n")
results %>%
  filter(outcome_label == "V1: supplemental",
         model == "primary+IPW") %>%
  select(exposure, term, rr_ci, p_fmt) %>%
  print(n = Inf)

cat("\n======================================================================\n")
cat("COMPARISON: Primary (no site) vs Sensitivity (with site)\n")
cat("======================================================================\n")
bind_rows(
  results %>%
    filter(outcome_label == "V1: primary (no site)",
           model == "primary+IPW") %>%
    mutate(model = "Primary (no site)"),
  results %>%
    filter(outcome_label == "V1: sensitivity (with site)",
           model == "site+IPW") %>%
    mutate(model = "With site")
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