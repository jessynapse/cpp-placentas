# ==============================================================================
# 01-load-clean.R
# CPP Placental Pathology → Infantile Hemangioma
# Updated June 2026
#
# Decisions:
#   - No GA restriction (Alexa twin paper approach)
#   - Twins excluded pending CHILDID merge (Alexa returns ~June 2026)
#   - MVM dedup via slice_max: 31 singleton conflicts checked,
#     only 3 in analytic sample, 0 IH cases — defensible
#   - Single IPW model only (LTFU non-differential by exposure, per Ellen)
#   - Chronic HTN from cpp_htn_20231031.sas7bdat
#   - Preterm, SGA, labor onset added per Ellen Francis
#   - SGA sex-specific per Freedman et al. (grouped by gest_age + infant_sex)
#   - Parity: na_tag() used to distinguish .P (primigravida → 0) from
#     .U (unknown → NA), per Chakraborti guidance (November 2023)
# ==============================================================================

library(haven)
library(tidyverse)
library(labelled)

setwd("/Users/wongjj/Downloads")

clean_numeric <- function(x) {
  x <- as.character(x)
  x[x %in% c(".u", ".U", "u", "U")] <- NA
  as.numeric(x)
}

# ------------------------------------------------------------------------------
# LOAD
# ------------------------------------------------------------------------------
path1    <- read_sas("path1all.sas7bdat")
path2    <- read_sas("path2all.sas7bdat")
ncppbasa <- read_sas("ncppbasa.sas7bdat")

diabetes_data <- read_sas("cpp_diabetes_20231026.sas7bdat") %>%
  select(MOMID, PREGID, DM)

htn_data <- read_sas("cpp_htn_20231031.sas7bdat") %>%
  select(MOMID, PREGID, HTN_PREG) %>%
  mutate(chronic_htn = as.integer(HTN_PREG == 1))

mvm_data <- read_csv("CPP_dataset_021226.csv", show_col_types = FALSE) %>%
  rename(PREGID = NINDB_NUMBER) %>%
  group_by(PREGID) %>%
  slice_max(MVM_SCORE, n = 1, with_ties = FALSE) %>%
  ungroup()

# ------------------------------------------------------------------------------
# MERGE
# ------------------------------------------------------------------------------
ncpp_vars <- ncppbasa %>% select(
  MOMID, PREGID, CHILDID,
  C1095, C1101,          # birthweight, gest age
  C983,                  # IH (no PWS)
  C31, C83, C81,         # age, weight, height
  C52, C50, C10,         # smoking, parity, plurality
  C275, C272,            # income, education
  C36, C554,             # marital, infant sex
  C303, C1092,           # race, vital status
  C1317,                 # labor onset
  SITE)

analytic <- path1 %>%
  inner_join(path2,         by = c("MOMID", "PREGID", "CHILDID")) %>%
  inner_join(ncpp_vars,     by = c("MOMID", "PREGID", "CHILDID")) %>%
  left_join(diabetes_data,  by = c("MOMID", "PREGID")) %>%
  left_join(htn_data,       by = c("MOMID", "PREGID")) %>%
  left_join(mvm_data,       by = "PREGID") %>%
  rename(birthweight = C1095, gest_age = C1101)

# ------------------------------------------------------------------------------
# OUTCOME + EXCLUSION
# ------------------------------------------------------------------------------
analytic <- analytic %>%
  mutate(
    IH_raw    = clean_numeric(C983),
    exclusion = ifelse(C1092 == 0 | (C1092 >= 51 & C1092 <= 58), 0, 1),
    ih1       = case_when(exclusion == 1 ~ NA_real_,
                          IH_raw == 2    ~ 0,
                          TRUE           ~ IH_raw),
    ih2       = case_when(exclusion == 1 ~ NA_real_,
                          IH_raw == 2    ~ 1,
                          TRUE           ~ IH_raw)
  )

# ------------------------------------------------------------------------------
# EXCLUDE TWINS + MISSING MVM
# ------------------------------------------------------------------------------
twin_pregids <- analytic %>%
  group_by(PREGID) %>% filter(n() > 1) %>%
  pull(PREGID) %>% unique()

analytic <- analytic %>%
  filter(!PREGID %in% twin_pregids, !is.na(MVM_SCORE))

# ------------------------------------------------------------------------------
# COVARIATES
# ------------------------------------------------------------------------------
analytic <- analytic %>%
  mutate(
    age     = clean_numeric(C31),
    bmi     = ifelse(clean_numeric(C81) * 0.0254 > 0,
                     (clean_numeric(C83) * 0.453592) /
                       (clean_numeric(C81) * 0.0254)^2, NA),
    smoking = case_when(
      clean_numeric(C52) == 0                             ~ 0,
      clean_numeric(C52) >= 1 & clean_numeric(C52) <= 19  ~ 1,
      clean_numeric(C52) >= 20                            ~ 2),
    
    # ------------------------------------------------------------------
    # PARITY — use na_tag() to distinguish .P (primigravida → 0)
    # from .U (unknown → NA), per Chakraborti guidance (November 2023)
    # na_tag() returns "p" for SAS .P and "u" for SAS .U
    # ------------------------------------------------------------------
    parity  = pmin(case_when(
      na_tag(C50) == "p" ~ 0,
      na_tag(C50) == "u" ~ NA_real_,
      is.na(C50)         ~ NA_real_,
      TRUE               ~ as.numeric(C50)), 4),
    
    plurality = case_when(
      clean_numeric(C10) == 1           ~ 1,
      clean_numeric(C10) %in% c(2,3,4) ~ 2),
    income  = clean_numeric(C275),
    educ    = case_when(
      clean_numeric(C272) %in% 0:2 ~ 0,
      clean_numeric(C272) == 3     ~ 1,
      clean_numeric(C272) == 4     ~ 2,
      clean_numeric(C272) %in% 5:6 ~ 3),
    marital = case_when(
      clean_numeric(C36) %in% c(2,3)   ~ 0,
      clean_numeric(C36) == 1          ~ 1,
      clean_numeric(C36) %in% c(4,5,6) ~ 2),
    infant_sex = case_when(
      clean_numeric(C554) == 1 ~ 1,
      clean_numeric(C554) == 2 ~ 0),
    race = case_when(
      clean_numeric(C303) == 3 ~ 8,
      TRUE                     ~ clean_numeric(C303)),
    site         = as.factor(SITE),
    dm           = DM,
    mvm          = MVM_SCORE,
    mvm3         = MVM_SCORE_CAT,
    mvm2         = as.integer(MVM_SCORE_CAT > 0),
    mvm_villous  = MVM_VIL,
    mvm_vascular = MVM_VES,
    ai           = AI_DI,
    ai3          = AI_3cat,
    chorangioma  = as.integer(PA201_56 == 1),
    preterm      = as.integer(gest_age < 37),
    preterm3     = case_when(
      gest_age < 32  ~ 2,
      gest_age < 37  ~ 1,
      gest_age >= 37 ~ 0),
    labor_onset  = case_when(
      clean_numeric(C1317) == 0 ~ "None",
      clean_numeric(C1317) == 1 ~ "Spontaneous",
      clean_numeric(C1317) == 2 ~ "Induced")
  )

# ------------------------------------------------------------------------------
# PARITY CHECK
# ------------------------------------------------------------------------------
cat("Parity distribution (should have low missingness ~1%):\n")
analytic %>%
  summarise(
    n_total   = n(),
    n_missing = sum(is.na(parity)),
    pct_miss  = round(100 * mean(is.na(parity)), 1),
    n_zero    = sum(parity == 0, na.rm = TRUE),
    pct_zero  = round(100 * mean(parity == 0, na.rm = TRUE), 1)
  ) %>%
  print()

# ------------------------------------------------------------------------------
# SGA — sex-specific per Freedman et al.
# ------------------------------------------------------------------------------
analytic <- analytic %>%
  group_by(gest_age, infant_sex) %>%
  mutate(
    bw_for_ga_z = (birthweight - mean(birthweight, na.rm = TRUE)) /
      sd(birthweight, na.rm = TRUE),
    SGA         = as.integer(bw_for_ga_z < -1.28)
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# QUICK CHECK — SGA by sex (should be ~10% in each group)
# ------------------------------------------------------------------------------
cat("SGA prevalence by sex (should be ~10% each):\n")
analytic %>%
  filter(!is.na(SGA), !is.na(infant_sex)) %>%
  group_by(infant_sex) %>%
  summarise(
    N       = n(),
    n_SGA   = sum(SGA == 1),
    pct_SGA = round(100 * mean(SGA == 1), 1),
    .groups = "drop"
  ) %>%
  print()

cat("\nSGA prevalence by MVM grade (high MVM should be ~38%):\n")
analytic %>%
  filter(!is.na(SGA), !is.na(mvm3)) %>%
  mutate(MVM = factor(mvm3, 0:2, c("None","Low","High"))) %>%
  group_by(MVM) %>%
  summarise(
    N       = n(),
    n_SGA   = sum(SGA == 1),
    pct_SGA = round(100 * mean(SGA == 1), 1),
    .groups = "drop"
  ) %>%
  print()

# ------------------------------------------------------------------------------
# SELECT + LABEL
# ------------------------------------------------------------------------------
analytic <- analytic %>%
  select(MOMID, PREGID, CHILDID,
         mvm, mvm3, mvm2, mvm_villous, mvm_vascular,
         ai, ai3, chorangioma,
         ih1, ih2, exclusion,
         age, bmi, smoking, parity, plurality,
         income, educ, marital, infant_sex, race, site,
         dm, chronic_htn,
         gest_age, birthweight, bw_for_ga_z, SGA,
         preterm, preterm3, labor_onset)

var_label(analytic) <- list(
  mvm          = "MVM continuous score (0-8)",
  mvm3         = "MVM grade: 0=None, 1=Low (1-2), 2=High (3+)",
  mvm2         = "MVM binary: 0=None, 1=Any",
  mvm_villous  = "MVM villous lesions: 0=No, 1=Yes",
  mvm_vascular = "MVM vascular lesions: 0=No, 1=Yes",
  ai           = "Acute inflammation: 0=No, 1=Yes",
  ai3          = "AI stage: 0=None, 1=Low, 2=High",
  chorangioma  = "Chorangioma (PA201_56) — pending Linda confirmation",
  ih1          = "IH at 1y V1 (primary): suspect=absent",
  ih2          = "IH at 1y V2 (sensitivity): suspect=present",
  exclusion    = "0=alive at 1y, 1=died or LTFU before 1y",
  age          = "Maternal age, years",
  bmi          = "Pre-pregnancy BMI, kg/m2",
  smoking      = "0=None, 1=1-19/day, 2=20+/day",
  parity       = "Prior viable births: 0/1/2/3/4+",
  plurality    = "1=Singleton, 2=Multiple",
  income       = "Family income (C275 raw)",
  educ         = "0=<HS, 1=Some HS, 2=HS grad, 3=College+",
  marital      = "0=Married, 1=Single, 2=Wid/Div/Sep",
  infant_sex   = "0=Female, 1=Male",
  race         = "1=White, 2=Black, 4=Puerto Rican, 8=Other",
  site         = "Study site (12 CPP sites)",
  dm           = "Pre-pregnancy diabetes: 0=No, 1=Yes",
  chronic_htn  = "Chronic hypertension (HTN_PREG==1): 0=No, 1=Yes",
  gest_age     = "Gestational age at delivery, weeks",
  birthweight  = "Birthweight, grams",
  bw_for_ga_z  = "Birthweight-for-GA z-score (within GA week x sex)",
  SGA          = "SGA <10th pctile for GA x sex: 0=No, 1=Yes",
  preterm      = "Preterm <37w: 0=No, 1=Yes",
  preterm3     = "0=Term, 1=Preterm 32-36w, 2=Very preterm <32w",
  labor_onset  = "Labor onset: None/Spontaneous/Induced")

# ------------------------------------------------------------------------------
# CHECK + SAVE
# ------------------------------------------------------------------------------
cat("\nN:", nrow(analytic), "| IH cases:", sum(analytic$ih1 == 1, na.rm=TRUE), "\n")

cat("\nCovariate missingness:\n")
for (v in c("age","bmi","smoking","parity","plurality","income",
            "educ","marital","infant_sex","race","dm","chronic_htn",
            "site","preterm","SGA")) {
  n_miss <- sum(is.na(analytic[[v]]))
  cat(sprintf("  %-15s %5d (%s%%)\n", v, n_miss,
              round(100*n_miss/nrow(analytic), 1)))
}

saveRDS(analytic, "analytic_sample.RDS")
write_csv(analytic, "analytic_sample.csv", na = "")
write_csv(tibble(variable = names(var_label(analytic)),
                 label    = unlist(var_label(analytic))), "labels.csv")

cat("\nSaved: analytic_sample.RDS / .csv / labels.csv\n")