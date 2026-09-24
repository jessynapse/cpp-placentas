# ==============================================================================
# 11-compare-original-corrected.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Side-by-side comparison of the original pipeline (01, 03, 04, 05) and the
# corrected pipeline (07, 08, 09, 10).
#
# Corrections in 07-10:
#   - 26 twins excluded (C10 == 1)
#   - Chronic HTN includes HTN_PREG level 4 (superimposed preeclampsia)
#   - gest_age and birthweight removed as imputation predictors, so no
#     covariates are left missing (original models silently dropped 36
#     infants with an observed outcome)
#   - Income reference = $4,000-$5,999 (does not change exposure estimates)
#
# Outputs (outputs/):
#   compare_sample.csv     sample sizes and key counts, original vs corrected
#   compare_estimates.csv  every estimate, original vs corrected
# ==============================================================================

library(tidyverse)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE)

fmt_n   <- function(x) format(x, big.mark = ",", trim = TRUE)
pct_n   <- function(x, n) sprintf("%.1f%% (%s/%s)", 100 * x / n, fmt_n(x), fmt_n(n))

# ------------------------------------------------------------------------------
# SAMPLE COMPARISON
# ------------------------------------------------------------------------------
orig <- readRDS(file.path(in_dir, "analytic_sample.RDS"))
corr <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS"))

# Infants actually used in adjusted models = observed outcome with complete
# covariates after imputation (imputation 1)
model_n <- function(file) {
  w <- readRDS(file.path(in_dir, file)) %>% filter(.imp == 1)
  covs <- c("age","bmi","race","educ","income","marital","smoking",
            "parity","dm","chronic_htn","infant_sex")
  obs <- w %>% filter(!is.na(ih1), exclusion == 0)
  sum(complete.cases(obs[, covs]))
}

sample_row <- function(d, n_model) {
  obs <- d %>% filter(!is.na(ih1), exclusion == 0)
  tibble(
    pregnancies            = fmt_n(nrow(d)),
    observed_ih_outcome    = pct_n(nrow(obs), nrow(d)),
    used_in_adjusted_models = fmt_n(n_model),
    definite_ih            = pct_n(sum(obs$ih1 == 1), nrow(obs)),
    ih_incl_suspect        = pct_n(sum(obs$ih2 == 1), nrow(obs)),
    mvm_any                = pct_n(sum(d$mvm2 == 1, na.rm = TRUE), sum(!is.na(d$mvm2))),
    ai_any                 = pct_n(sum(d$ai == 1, na.rm = TRUE), sum(!is.na(d$ai))),
    chronic_htn            = pct_n(sum(d$chronic_htn == 1, na.rm = TRUE),
                                   sum(!is.na(d$chronic_htn))),
    died_before_1y         = pct_n(sum(d$exclusion == 1), nrow(d)),
    lost_to_follow_up      = pct_n(sum(d$exclusion == 0 & is.na(d$ih1)), nrow(d))
  )
}

compare_sample <- bind_rows(
  sample_row(orig, model_n("weighted_long.RDS"))           %>% mutate(pipeline = "Original (01-05)"),
  sample_row(corr, model_n("weighted_long_corrected.RDS")) %>% mutate(pipeline = "Corrected (07-10)")
) %>%
  pivot_longer(-pipeline, names_to = "measure") %>%
  pivot_wider(names_from = pipeline, values_from = value)

# ------------------------------------------------------------------------------
# ESTIMATE COMPARISON
# Original labels from 05 mapped to the 10 labels
# ------------------------------------------------------------------------------
r_orig <- readRDS(file.path(in_dir, "outcome_results.RDS")) %>%
  mutate(
    outcome = ifelse(grepl("^V1", outcome_label), "ih1", "ih2"),
    analysis_type = case_when(
      grepl("crude", outcome_label)        ~ "crude",
      grepl("supplemental", outcome_label) ~ "supplemental",
      grepl("with site", outcome_label)    ~ "site",
      TRUE                                 ~ "primary"),
    model = recode(model,
                   "primary"       = "adjusted",
                   "primary+IPW"   = "adjusted+IPW",
                   "site adjusted" = "adjusted+site",
                   "site+IPW"      = "adjusted+site+IPW"),
    exposure = recode(exposure, "MVM x AI joint" = "MVM x AI interaction")
  ) %>%
  select(outcome, analysis_type, exposure, model, term, orig = rr_ci)

r_corr <- readRDS(file.path(in_dir, "outcome_results_corrected.RDS")) %>%
  mutate(analysis_type = str_extract(analysis, "crude|supplemental|site|primary")) %>%
  select(outcome, analysis_type, exposure, model, term,
         n_cases, n_obs, corrected = rr_ci, rr)

compare_est <- full_join(r_orig, r_corr,
                         by = c("outcome", "analysis_type", "exposure", "model", "term")) %>%
  mutate(outcome = recode(outcome, ih1 = "V1 definite", ih2 = "V2 suspect = present"),
         orig = replace_na(orig, "not fitted in 05")) %>%
  arrange(outcome, factor(analysis_type, c("crude","primary","site","supplemental")),
          exposure, model, term) %>%
  select(outcome, analysis_type, exposure, model, term,
         original = orig, corrected, n_cases, n_obs)

# ------------------------------------------------------------------------------
# PRINT + SAVE
# ------------------------------------------------------------------------------
cat("\n==================== SAMPLE ====================\n")
print(compare_sample, n = Inf, width = Inf)

cat("\n==================== KEY ESTIMATES (adjusted + IPW) ====================\n")
compare_est %>%
  filter(model %in% c("crude", "adjusted+IPW", "adjusted+site+IPW"),
         analysis_type != "supplemental") %>%
  print(n = Inf, width = Inf)

cat("\n==================== SUPPLEMENTAL (adjusted + IPW) ====================\n")
compare_est %>%
  filter(analysis_type == "supplemental", model == "adjusted+IPW") %>%
  print(n = Inf, width = Inf)

write_csv(compare_sample, file.path(out_dir, "compare_sample.csv"))
write_csv(compare_est,    file.path(out_dir, "compare_estimates.csv"))
cat("\nSaved: outputs/compare_sample.csv, outputs/compare_estimates.csv\n")
