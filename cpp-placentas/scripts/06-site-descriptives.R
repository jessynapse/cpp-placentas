# ==============================================================================
# 06-site-descriptives.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Descriptive analysis by study site (12 CPP centres), prompted by:
#   - Alexa Freedman comment: pathology varies by site more than expected,
#     does IH prevalence / ascertainment also vary by site?
#   - Script 05 site-adjusted sensitivity model attenuating all three
#     primary estimates toward the null
#
# Uses analytic_sample.RDS from 01-load-clean.R (current sample, N = 45,296,
# which still includes 26 twins and chronic_htn coded as HTN_PREG == 1 only)
#
# Site codes are numeric. The SAS format catalogue with site names is not
# available. Site 5 = Boston (per 05-outcomemodels.R).
#
# Outputs (outputs/):
#   site_A_sample_followup.csv   N, deaths, LTFU, observed outcome by site
#   site_B_outcome.csv           definite and suspect IH by site
#   site_C_exposure.csv          MVM and AI prevalence by site
#   site_D_population.csv        race and maternal characteristics by site
#   site_E_ih_by_exposure.csv    IH by MVM and by AI within each site,
#                                with crude within-site RRs
#
# No p-values. Percentages reported before counts.
# ==============================================================================

library(tidyverse)

in_dir  <- "/Users/wongjj/Downloads"   # folder with analytic_sample.RDS
out_dir <- "outputs"                   # relative to the project folder
dir.create(out_dir, showWarnings = FALSE)

df <- readRDS(file.path(in_dir, "analytic_sample.RDS")) %>%
  mutate(
    site     = factor(as.character(site),
                      levels = as.character(sort(unique(as.numeric(as.character(site)))))),
    ih_obs   = !is.na(ih1) & exclusion == 0,
    died     = exclusion == 1,
    ltfu     = exclusion == 0 & is.na(ih1)
  )

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------
# "2.0% (225/11,334)"
pct_n <- function(x, n, digits = 1) {
  ifelse(n == 0, "NA",
         sprintf(paste0("%.", digits, "f%% (%s/%s)"),
                 100 * x / n, format(x, big.mark = ","), format(n, big.mark = ",")))
}

# Crude RR with Wald 95% CI on the log scale, formatted "1.25 (1.07, 1.46)"
crude_rr <- function(a, n1, c, n0) {
  if (a == 0 || c == 0) return("not estimable (0 cases in a group)")
  rr <- (a / n1) / (c / n0)
  se <- sqrt(1/a - 1/n1 + 1/c - 1/n0)
  sprintf("%.2f (%.2f, %.2f)", rr, exp(log(rr) - 1.96 * se), exp(log(rr) + 1.96 * se))
}

add_total <- function(data, fun) {
  bind_rows(fun(data) %>% mutate(site = as.character(site)),
            fun(data %>% mutate(site = "All sites")))
}

# ------------------------------------------------------------------------------
# A. SAMPLE SIZE AND FOLLOW-UP BY SITE
# ------------------------------------------------------------------------------
tab_A <- add_total(df, function(d) d %>%
  group_by(site) %>%
  summarise(
    n_pregnancies   = n(),
    pct_of_sample   = sprintf("%.1f%%", 100 * n() / nrow(df)),
    died_before_1y  = pct_n(sum(died), n()),
    lost_to_fu      = pct_n(sum(ltfu), n()),
    observed_ih     = pct_n(sum(ih_obs), n()),
    .groups = "drop"))

# ------------------------------------------------------------------------------
# B. OUTCOME BY SITE (among infants with an observed IH outcome)
# ------------------------------------------------------------------------------
tab_B <- add_total(df, function(d) d %>%
  filter(ih_obs) %>%
  group_by(site) %>%
  summarise(
    n_observed        = n(),
    ih_definite       = pct_n(sum(ih1 == 1), n()),
    ih_suspect_only   = pct_n(sum(ih2 == 1 & ih1 == 0), n()),
    ih_incl_suspect   = pct_n(sum(ih2 == 1), n()),
    share_of_all_definite_cases =
      sprintf("%.1f%%", 100 * sum(ih1 == 1) / sum(df$ih1 == 1, na.rm = TRUE)),
    .groups = "drop"))

# ------------------------------------------------------------------------------
# C. EXPOSURE PREVALENCE BY SITE (all pregnancies)
# ------------------------------------------------------------------------------
tab_C <- add_total(df, function(d) d %>%
  group_by(site) %>%
  summarise(
    n_pregnancies  = n(),
    mvm_any        = pct_n(sum(mvm2 == 1, na.rm = TRUE), sum(!is.na(mvm2))),
    mvm_high_grade = pct_n(sum(mvm3 == 2, na.rm = TRUE), sum(!is.na(mvm3))),
    mvm_score_mean = sprintf("%.2f (%.2f)", mean(mvm, na.rm = TRUE), sd(mvm, na.rm = TRUE)),
    ai_any         = pct_n(sum(ai == 1, na.rm = TRUE), sum(!is.na(ai))),
    ai_low_stage   = pct_n(sum(ai3 == 1, na.rm = TRUE), sum(!is.na(ai3))),
    ai_high_stage  = pct_n(sum(ai3 == 2, na.rm = TRUE), sum(!is.na(ai3))),
    .groups = "drop"))

# ------------------------------------------------------------------------------
# D. POPULATION CHARACTERISTICS BY SITE
# ------------------------------------------------------------------------------
tab_D <- add_total(df, function(d) d %>%
  group_by(site) %>%
  summarise(
    n_pregnancies  = n(),
    white          = pct_n(sum(race == 1, na.rm = TRUE), sum(!is.na(race))),
    black          = pct_n(sum(race == 2, na.rm = TRUE), sum(!is.na(race))),
    puerto_rican   = pct_n(sum(race == 4, na.rm = TRUE), sum(!is.na(race))),
    other          = pct_n(sum(race == 8, na.rm = TRUE), sum(!is.na(race))),
    maternal_age   = sprintf("%.1f (%.1f)", mean(age, na.rm = TRUE), sd(age, na.rm = TRUE)),
    preterm        = pct_n(sum(preterm == 1, na.rm = TRUE), sum(!is.na(preterm))),
    .groups = "drop"))

# ------------------------------------------------------------------------------
# E. IH BY EXPOSURE WITHIN EACH SITE (observed outcome only)
# Crude within-site RRs are descriptive only. Small cells make many unstable.
# ------------------------------------------------------------------------------
within_site <- function(d, exp_var) {
  d %>%
    filter(ih_obs, !is.na(.data[[exp_var]])) %>%
    group_by(site) %>%
    summarise(
      a  = sum(ih1 == 1 & .data[[exp_var]] == 1), n1 = sum(.data[[exp_var]] == 1),
      c  = sum(ih1 == 1 & .data[[exp_var]] == 0), n0 = sum(.data[[exp_var]] == 0),
      .groups = "drop") %>%
    rowwise() %>%
    mutate(ih_exposed   = pct_n(a, n1),
           ih_unexposed = pct_n(c, n0),
           crude_rr     = crude_rr(a, n1, c, n0)) %>%
    ungroup()
}

tab_E <- bind_rows(
  add_total(df, function(d) within_site(d, "mvm2")) %>% mutate(exposure = "MVM any vs none"),
  add_total(df, function(d) within_site(d, "ai"))   %>% mutate(exposure = "AI any vs none")
) %>%
  select(exposure, site, n_cases = a, ih_exposed, ih_unexposed, crude_rr, n1, c, n0) %>%
  mutate(n_cases = n_cases + c, n_observed = n1 + n0) %>%
  select(exposure, site, n_cases, n_observed, ih_exposed, ih_unexposed, crude_rr)

# ------------------------------------------------------------------------------
# PRINT + SAVE
# ------------------------------------------------------------------------------
for (nm in c("A", "B", "C", "D", "E")) {
  cat("\n==============================================================\n")
  cat("TABLE", nm, "\n")
  cat("==============================================================\n")
  print(get(paste0("tab_", nm)), n = Inf, width = Inf)
}

write_csv(tab_A, file.path(out_dir, "site_A_sample_followup.csv"))
write_csv(tab_B, file.path(out_dir, "site_B_outcome.csv"))
write_csv(tab_C, file.path(out_dir, "site_C_exposure.csv"))
write_csv(tab_D, file.path(out_dir, "site_D_population.csv"))
write_csv(tab_E, file.path(out_dir, "site_E_ih_by_exposure.csv"))

cat("\nSaved 5 tables to", out_dir, "\n")
