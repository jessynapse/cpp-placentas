# ==============================================================================
# 00_pipeline.R
# CPP PLACENTAL PATHOLOGY → INFANTILE HEMANGIOMA
# Updated May 2026
#
# Key decisions:
#   - No GA restriction (consistent with Alexa's twin paper)
#   - Twins excluded pending CHILDID merge with Alexa (~June 2026)
#   - MVM dedup: slice_max (take higher score) — NOT removing conflicts
#     Reason: 31 singleton conflicts checked, only 3 in analytic sample,
#     0 IH cases — keeping them with higher score is defensible
#   - Single IPW model only (LTFU non-differential by exposure per Ellen)
#   - Chronic HTN added (missing from this pipeline)
#   - Preterm birth added (Ellen's suggestion)
#   - SGA/birthweight for GA added (Ellen's suggestion)
#   - C1317 labor onset added (spontaneous vs induced)
#   - HTN loaded from cpp_htn_20231031.sas7bdat
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
path1         <- read_sas("path1all.sas7bdat")
path2         <- read_sas("path2all.sas7bdat")
ncppbasa      <- read_sas("ncppbasa.sas7bdat")
diabetes_data <- read_sas("cpp_diabetes_20231026.sas7bdat") %>%
  select(MOMID, PREGID, DM)
htn_data      <- read_sas("cpp_htn_20231031.sas7bdat") %>%
  select(MOMID, PREGID, HTN_PREG) %>%
  mutate(chronic_htn = as.integer(HTN_PREG == 1))
mvm_data      <- read_csv("CPP_dataset_021226.csv", show_col_types = FALSE) %>%
  rename(PREGID = NINDB_NUMBER)

# ------------------------------------------------------------------------------
# DEDUPLICATE MVM
# Take higher MVM_SCORE row per PREGID (twins + 31 singleton conflicts)
# 31 singleton conflicts: only 3 in analytic sample, 0 IH cases — kept
# Pending Alexa (~June 2026): confirm approach for singleton conflicts
# ------------------------------------------------------------------------------
mvm_clean <- mvm_data %>%
  group_by(PREGID) %>%
  slice_max(MVM_SCORE, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("MVM dedup:", nrow(mvm_data), "->", nrow(mvm_clean), "rows\n")

# ------------------------------------------------------------------------------
# MERGE — no GA restriction
# ------------------------------------------------------------------------------
dat <- path1 %>%
  inner_join(path2, by = c("MOMID", "PREGID", "CHILDID")) %>%
  inner_join(
    ncppbasa %>% select(
      MOMID, PREGID, CHILDID,
      C1095, C1101,   # birthweight, gestational age
      C983,           # IH outcome (no PWS)
      C31, C83, C81,  # age, weight, height
      C52, C50, C10,  # smoking, parity, plurality
      C275, C272,     # income, education
      C36, C554,      # marital, infant sex
      C303, C1092,    # race, vital status
      C1317,          # labor onset: 1=spontaneous, 2=induced
      SITE
    ),
    by = c("MOMID", "PREGID", "CHILDID")
  ) %>%
  left_join(diabetes_data, by = c("MOMID", "PREGID")) %>%
  left_join(htn_data,      by = c("MOMID", "PREGID")) %>%
  left_join(mvm_clean,     by = "PREGID") %>%
  rename(
    birthweight = C1095,
    gest_age    = C1101
  )

cat("After merge:", nrow(dat), "rows,", ncol(dat), "cols\n")

# ------------------------------------------------------------------------------
# OUTCOME AND EXCLUSION
# ------------------------------------------------------------------------------
dat <- dat %>%
  mutate(
    IH_raw    = clean_numeric(C983),
    mat_age   = clean_numeric(C31),
    exclusion = ifelse(C1092 == 0 | (C1092 >= 51 & C1092 <= 58), 0, 1),
    ih1 = case_when(      # primary: suspect -> absent
      exclusion == 1 ~ NA_real_,
      IH_raw == 2    ~ 0,
      TRUE           ~ IH_raw
    ),
    ih2 = case_when(      # sensitivity: suspect -> present
      exclusion == 1 ~ NA_real_,
      IH_raw == 2    ~ 1,
      TRUE           ~ IH_raw
    )
  )

# ------------------------------------------------------------------------------
# EXCLUDE TWINS
# ------------------------------------------------------------------------------
twin_pregids <- dat %>%
  group_by(PREGID) %>%
  filter(n() > 1) %>%
  pull(PREGID) %>%
  unique()

dat_singletons <- dat %>%
  filter(!PREGID %in% twin_pregids)

cat("Twins excluded:", length(twin_pregids), "PREGIDs,",
    sum(dat$PREGID %in% twin_pregids), "rows\n")
cat("After twin exclusion:", nrow(dat_singletons), "\n")

# ------------------------------------------------------------------------------
# ANALYTIC SAMPLE
# ------------------------------------------------------------------------------
analytic <- dat_singletons %>%
  filter(!is.na(MVM_SCORE))

cat("======================================================================\n")
cat("ANALYTIC SAMPLE\n")
cat("======================================================================\n")
cat("N:                  ", nrow(analytic), "\n")
cat("Non-survivors:      ", sum(analytic$exclusion == 1, na.rm=TRUE), "\n")
cat("Survivors LTFU:     ",
    sum(analytic$exclusion == 0 & is.na(analytic$ih1)), "\n")
cat("IH observed:        ",
    sum(analytic$exclusion == 0 & !is.na(analytic$ih1)), "\n")
cat("IH cases (V1):      ", sum(analytic$ih1 == 1, na.rm=TRUE), "\n")
cat("IH prevalence:      ",
    round(100 * mean(analytic$ih1 == 1, na.rm=TRUE), 2), "%\n")

# ------------------------------------------------------------------------------
# COVARIATES
# ------------------------------------------------------------------------------
analytic <- analytic %>%
  mutate(
    age           = mat_age,
    weight_kg     = clean_numeric(C83) * 0.453592,
    height_m      = clean_numeric(C81) * 0.0254,
    bmi           = ifelse(height_m > 0, weight_kg / height_m^2, NA),
    
    smoking_raw   = clean_numeric(C52),
    smoking       = case_when(
      is.na(smoking_raw)                   ~ NA_real_,
      smoking_raw == 0                     ~ 0,
      smoking_raw >= 1 & smoking_raw <= 19 ~ 1,
      smoking_raw >= 20                    ~ 2
    ),
    
    parity_raw    = case_when(
      as.character(C50) == "p" ~ 0,
      as.character(C50) == "u" ~ NA_real_,
      TRUE ~ as.numeric(as.character(C50))
    ),
    parity        = ifelse(parity_raw >= 4, 4, parity_raw),
    
    plurality     = case_when(
      clean_numeric(C10) == 1            ~ 1,
      clean_numeric(C10) %in% c(2, 3, 4) ~ 2,
      TRUE                               ~ NA_real_
    ),
    
    income        = clean_numeric(C275),
    
    # Education: 0=<HS, 1=Some HS (ref), 2=HS grad, 3=Some college+
    # Note: Some HS is reference category (most common in CPP)
    educ          = case_when(
      clean_numeric(C272) %in% 0:2 ~ 0,
      clean_numeric(C272) == 3      ~ 1,
      clean_numeric(C272) == 4      ~ 2,
      clean_numeric(C272) %in% 5:6  ~ 3,
      TRUE                          ~ NA_real_
    ),
    
    marital       = case_when(
      clean_numeric(C36) %in% c(2, 3)    ~ 0,  # married/common law (ref)
      clean_numeric(C36) == 1            ~ 1,  # single
      clean_numeric(C36) %in% c(4, 5, 6) ~ 2,  # widowed/divorced/sep
      TRUE                               ~ NA_real_
    ),
    
    infant_sex    = case_when(
      clean_numeric(C554) == 1 ~ 1,   # male
      clean_numeric(C554) == 2 ~ 0,   # female (ref)
      TRUE                     ~ NA_real_
    ),
    
    race          = case_when(
      clean_numeric(C303) == 3 ~ 8,
      TRUE                     ~ clean_numeric(C303)
    ),
    
    site          = as.factor(SITE),
    dm            = DM,
    # chronic_htn already constructed from htn_data join
    
    # Exposures
    mvm           = MVM_SCORE,
    mvm3          = MVM_SCORE_CAT,
    mvm2          = as.integer(MVM_SCORE_CAT > 0),
    mvm_villous   = MVM_VIL,
    mvm_vascular  = MVM_VES,
    ai            = AI_DI,
    ai3           = AI_3cat,
    
    chorangioma   = as.integer(PA201_56 == 1),
    
    # Perinatal outcomes (Ellen's suggestions)
    preterm       = as.integer(gest_age < 37),
    preterm3      = case_when(
      gest_age < 32  ~ 2,   # very preterm
      gest_age < 37  ~ 1,   # preterm
      gest_age >= 37 ~ 0,   # term
      TRUE           ~ NA_real_
    ),
    labor_onset   = case_when(
      clean_numeric(C1317) == 1 ~ "Spontaneous",
      clean_numeric(C1317) == 2 ~ "Induced",
      clean_numeric(C1317) == 0 ~ "None",
      TRUE                      ~ NA_character_
    )
  ) %>%
  # Birthweight for gestational age z-score (within GA week)
  group_by(gest_age) %>%
  mutate(
    bw_for_ga_z = (birthweight - mean(birthweight, na.rm=TRUE)) /
      sd(birthweight, na.rm=TRUE),
    SGA         = as.integer(bw_for_ga_z < -1.28)
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# SELECT FINAL VARIABLES
# ------------------------------------------------------------------------------
analytic <- analytic %>%
  select(
    MOMID, PREGID, CHILDID,
    # Exposures
    mvm, mvm3, mvm2, mvm_villous, mvm_vascular,
    ai, ai3,
    chorangioma,
    # Outcomes
    ih1, ih2,
    exclusion,
    # Covariates
    age, bmi,
    smoking, parity, plurality,
    income, educ, marital,
    infant_sex, race, site,
    dm, chronic_htn,
    # Perinatal (Ellen)
    gest_age, birthweight,
    bw_for_ga_z, SGA,
    preterm, preterm3, labor_onset
  )

# ------------------------------------------------------------------------------
# VARIABLE LABELS
# ------------------------------------------------------------------------------
var_label(analytic) <- list(
  mvm           = "MVM continuous score (0-8)",
  mvm3          = "MVM grade: 0=None, 1=Low (1-2), 2=High (3+)",
  mvm2          = "MVM binary: 0=None, 1=Any",
  mvm_villous   = "MVM villous lesions: 0=No, 1=Yes",
  mvm_vascular  = "MVM vascular lesions: 0=No, 1=Yes",
  ai            = "Acute inflammation: 0=No, 1=Yes",
  ai3           = "AI stage: 0=None, 1=Low, 2=High",
  chorangioma   = "Chorangioma (PA201_56): 0=No, 1=Yes — pending Linda confirmation",
  ih1           = "IH at 1y V1 (primary): suspect=absent",
  ih2           = "IH at 1y V2 (sensitivity): suspect=present",
  exclusion     = "0=alive/observed at 1y, 1=died or LTFU before 1y",
  age           = "Maternal age, years",
  bmi           = "Pre-pregnancy BMI, kg/m2",
  smoking       = "Smoking: 0=None, 1=1-19/day, 2=20+/day",
  parity        = "Prior viable births: 0/1/2/3/4+",
  plurality     = "1=Singleton, 2=Multiple",
  income        = "Family income category (C275)",
  educ          = "Education: 0=<HS, 1=Some HS, 2=HS grad, 3=College+",
  marital       = "Marital: 0=Married, 1=Single, 2=Wid/Div/Sep",
  infant_sex    = "0=Female, 1=Male",
  race          = "1=White, 2=Black, 4=Puerto Rican, 8=Other",
  site          = "Study site (12 CPP sites)",
  dm            = "Pre-pregnancy diabetes: 0=No, 1=Yes",
  chronic_htn   = "Chronic hypertension (HTN_PREG==1): 0=No, 1=Yes",
  gest_age      = "Gestational age at delivery, weeks",
  birthweight   = "Birthweight, grams",
  bw_for_ga_z   = "Birthweight-for-GA z-score (within GA week)",
  SGA           = "Small for gestational age (<10th pctile): 0=No, 1=Yes",
  preterm       = "Preterm birth (<37w): 0=No, 1=Yes",
  preterm3      = "0=Term, 1=Preterm 32-36w, 2=Very preterm <32w",
  labor_onset   = "Labor onset: Spontaneous / Induced / None"
)

# ------------------------------------------------------------------------------
# MISSINGNESS CHECK
# ------------------------------------------------------------------------------
cat("\nCovariate missingness:\n")
covs <- c("age", "bmi", "smoking", "parity", "plurality", "income",
          "educ", "marital", "infant_sex", "race", "dm",
          "chronic_htn", "site", "preterm", "SGA")
for (v in covs) {
  n_miss <- sum(is.na(analytic[[v]]))
  pct    <- round(100 * n_miss / nrow(analytic), 1)
  cat(sprintf("  %-15s %5d missing (%s%%)\n", v, n_miss, pct))
}

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
saveRDS(analytic, "analytic_sample.RDS")
write_csv(analytic, "analytic_sample.csv", na = "")

labels_df <- tibble(
  variable = names(var_label(analytic)),
  label    = unlist(var_label(analytic))
)
write_csv(labels_df, "labels.csv")

cat("\nSaved: analytic_sample.RDS, analytic_sample.csv, labels.csv\n")
cat("N:", nrow(analytic), "| IH cases:", sum(analytic$ih1 == 1, na.rm=TRUE), "\n")