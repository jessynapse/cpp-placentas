# ==============================================================================
# 16-joint-mvm-ai.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# MVM and AI together, requested by Alexa Freedman ("what happens if you do
# put both in the model?") and Linda Ernst ("the effects of MVM and AI
# together"). Needed before the paper can call the associations
# "independent".
#
#   1. Mutually adjusted: mvm2 + ai in the same model
#   2. Four-level joint exposure: neither (reference), MVM only, AI only,
#      both
#
# Corrected pipeline (07-10): 40,700 infants, 700 definite IH cases,
# 20 imputations, Rubin's rules, modified Poisson GEE (exchangeable,
# clustered on MOMID). Weights from 09 (denominator includes mvm2 and ai).
# Models: crude, secondary (11 covariates + IPW), primary (+ site + IPW).
# Definite IH only.
#
# Outputs (outputs/):
#   joint_models.csv    estimates with case and observation counts
#   joint_counts.csv    IH by joint group, and MVM-AI overlap
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

fmt <- function(x, n) sprintf("%.1f%% (%s/%s)", 100 * x / n,
                              format(x, big.mark = ",", trim = TRUE),
                              format(n, big.mark = ",", trim = TRUE))

# ------------------------------------------------------------------------------
# DATA (as 10), trimmed to model columns to limit memory
# ------------------------------------------------------------------------------
weighted <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS"))

# Overlap in the full analytic sample (45,268 pregnancies, not imputed)
base <- weighted %>% filter(.imp == 0)
overlap <- tibble(
  measure = c("Any AI among pregnancies with any MVM",
              "Any MVM among pregnancies with any AI",
              "Any AI among pregnancies without MVM",
              "Both MVM and AI, all pregnancies"),
  value = c(fmt(sum(base$mvm2 == 1 & base$ai == 1), sum(base$mvm2 == 1)),
            fmt(sum(base$mvm2 == 1 & base$ai == 1), sum(base$ai == 1)),
            fmt(sum(base$mvm2 == 0 & base$ai == 1), sum(base$mvm2 == 0)),
            fmt(sum(base$mvm2 == 1 & base$ai == 1), nrow(base))))
rm(base)

imp_data <- weighted %>%
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
                                 mvm2 == 1 & ai == 1 ~ "Both MVM and AI"),
                       levels = c("Neither", "MVM only", "AI only", "Both MVM and AI")),
    mvm2      = factor(mvm2),
    ai        = factor(ai)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor)) %>%
  select(.imp, MOMID, ih1, ipw, mvm2, ai, joint, site, age, bmi, race, educ,
         income, marital, smoking, parity, dm, chronic_htn, infant_sex) %>%
  arrange(MOMID) %>%
  split(.$.imp)
rm(weighted); invisible(gc())

# ------------------------------------------------------------------------------
# MODELS
# ------------------------------------------------------------------------------
pool_gee <- function(fits) {
  m        <- length(fits)
  coef_mat <- do.call(rbind, lapply(fits, coef))
  se_mat   <- do.call(rbind, lapply(fits, function(f) sqrt(diag(vcov(f)))))
  Q_bar <- colMeans(coef_mat)
  SE    <- sqrt(colMeans(se_mat^2) + (1 + 1/m) * apply(coef_mat, 2, var))
  tibble(term = names(Q_bar), rr = exp(Q_bar),
         lci = exp(Q_bar - 1.96 * SE), uci = exp(Q_bar + 1.96 * SE), n_imps = m)
}

confounders <- "age + bmi + race + educ + income + marital +
                smoking + parity + dm + chronic_htn + infant_sex"
models <- c("Crude" = NA, "Secondary: covariates + IPW, no site" = "",
            "Primary: covariates + site + IPW" = "+ site")

run <- function(exposure, extra) {
  rhs <- if (is.na(extra)) exposure else paste(exposure, "+", confounders, extra)
  f <- as.formula(paste("ih1 ~", rhs))
  use_ipw <- !is.na(extra)
  fit_one <- function(d) {
    d$.w <- if (use_ipw) d$ipw else rep(1, nrow(d))
    tryCatch(geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                    corstr = "exchangeable", weights = .w),
             error = function(err) err)
  }
  fits <- mclapply(imp_data, fit_one, mc.cores = n_cores)
  redo <- which(!sapply(fits, function(x) inherits(x, "geeglm")))
  for (i in redo) fits[[i]] <- fit_one(imp_data[[i]])
  bad <- sapply(fits, function(x) !inherits(x, "geeglm"))
  if (any(bad)) message("  ", sum(bad), " fit(s) failed: ", rhs)
  pool_gee(fits[!bad]) %>% filter(grepl("^mvm2|^ai|^joint", term))
}

specs <- c("Mutually adjusted (MVM + AI)" = "mvm2 + ai",
           "Joint 4-level exposure"        = "joint")

results <- imap_dfr(specs, function(expo, spec) imap_dfr(models, function(extra, lab) {
  cat(" ", spec, "/", lab, "\n")
  run(expo, extra) %>% mutate(analysis = spec, model = lab)
})) %>%
  mutate(
    contrast = recode(term,
                      mvm21 = "MVM any vs none (adjusted for AI)",
                      ai1   = "AI any vs none (adjusted for MVM)",
                      `jointMVM only` = "MVM only vs neither",
                      `jointAI only`  = "AI only vs neither",
                      `jointBoth MVM and AI` = "Both vs neither"),
    n_cases = sum(imp_data[[1]]$ih1 == 1), n_obs = nrow(imp_data[[1]]),
    rr_ci   = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci)) %>%
  select(analysis, contrast, model, n_cases, n_obs, rr_ci, rr, lci, uci, n_imps)

# IH by joint group (outcome and exposures never imputed)
counts <- imp_data[[1]] %>%
  group_by(joint) %>%
  summarise(n_cases = sum(ih1 == 1), n_obs = n(), .groups = "drop") %>%
  mutate(ih_pct = fmt(n_cases, n_obs))

# ------------------------------------------------------------------------------
# PRINT + SAVE
# ------------------------------------------------------------------------------
cat("\n==================== MVM AND AI OVERLAP (45,268 pregnancies) ====================\n")
print(overlap, width = Inf)
cat("\n==================== IH BY JOINT GROUP ====================\n")
print(counts, width = Inf)
cat("\n==================== ESTIMATES (definite IH) ====================\n")
print(results %>% select(analysis, contrast, model, rr_ci), n = Inf, width = Inf)

write_csv(results, file.path(out_dir, "joint_models.csv"))
write_csv(bind_rows(overlap %>% mutate(table = "overlap"),
                    counts %>% transmute(table = "IH by joint group",
                                         measure = as.character(joint), value = ih_pct)),
          file.path(out_dir, "joint_counts.csv"))
cat("\nSaved: outputs/joint_models.csv, outputs/joint_counts.csv\n")
