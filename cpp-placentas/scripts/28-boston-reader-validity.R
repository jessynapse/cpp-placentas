# ==============================================================================
# 28-boston-reader-validity.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Stefanie Hinkle: "the question is if we think that they were reviewed
# correctly at Boston by this single person?"
#
# A single reader gives consistency, not necessarily accuracy. As a check of
# construct validity, compare expected associations of the placental
# lesions with birth outcomes in Boston (92.6% one pathologist) vs the other
# 11 sites. If Boston readings are valid, expected associations should be at
# least as strong in Boston:
#   MVM -> small for gestational age, preterm birth, preeclampsia
#   AI  -> preterm birth, spontaneous labor
# Crude risk ratios with Wald 95% CIs, all 45,268 pregnancies with
# non-missing values (outcome at age 1 not required).
# Preeclampsia = HTN_PREG 3 (preeclampsia) or 4 (superimposed).
#
# Output: outputs/boston_reader_validity.csv. No p-values.
# ==============================================================================

library(tidyverse)
library(haven)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")

h <- read_sas(file.path(in_dir, "cpp_htn_20231031.sas7bdat")) %>%
  zap_labels() %>% select(MOMID, PREGID, HTN_PREG)

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  left_join(h, by = c("MOMID", "PREGID")) %>%
  mutate(group = ifelse(as.character(site) == "5", "Boston", "Other 11 sites"),
         preeclampsia = as.integer(HTN_PREG %in% c(3, 4)),
         spontaneous_labor = as.integer(labor_onset == "Spontaneous"))

rr <- function(dd, out, expo) {
  dd <- dd[!is.na(dd[[out]]) & !is.na(dd[[expo]]), ]
  a <- sum(dd[[out]] == 1 & dd[[expo]] == 1); n1 <- sum(dd[[expo]] == 1)
  c <- sum(dd[[out]] == 1 & dd[[expo]] == 0); n0 <- sum(dd[[expo]] == 0)
  r <- (a / n1) / (c / n0); se <- sqrt(1/a - 1/n1 + 1/c - 1/n0)
  tibble(exposed = sprintf("%.1f%% (%s/%s)", 100 * a / n1, a, n1),
         unexposed = sprintf("%.1f%% (%s/%s)", 100 * c / n0, c, n0),
         rr_ci = sprintf("%.2f (%.2f, %.2f)", r, exp(log(r) - 1.96 * se), exp(log(r) + 1.96 * se)))
}

checks <- tribble(
  ~outcome,             ~exposure, ~label,
  "SGA",                "mvm2",    "Small for gestational age, by MVM",
  "preterm",            "mvm2",    "Preterm birth, by MVM",
  "preeclampsia",       "mvm2",    "Preeclampsia, by MVM",
  "preterm",            "ai",      "Preterm birth, by AI",
  "spontaneous_labor",  "ai",      "Spontaneous labor, by AI")

out <- pmap_dfr(checks, function(outcome, exposure, label)
  map_dfr(c("Boston", "Other 11 sites"), function(g)
    rr(d[d$group == g, ], outcome, exposure) %>% mutate(check = label, group = g))) %>%
  select(check, group, exposed, unexposed, rr_ci)

print(out, n = Inf, width = Inf)
write_csv(out, file.path(out_dir, "boston_reader_validity.csv"))
cat("\nSaved: outputs/boston_reader_validity.csv\n")
