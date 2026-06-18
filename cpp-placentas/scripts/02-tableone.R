# ==============================================================================
# 02-tableone.R
# CPP Placental Pathology → Infantile Hemangioma
# Updated June 2026
# Produces: Table 1 by MVM any/none and AI any/none (binary primary)
#           Supplemental table by IH status
#           Cross-tabs for results section
#           Differential exclusion/LTFU check
# ==============================================================================

library(tidyverse)
library(table1)
library(tableone)
library(labelled)

setwd("/Users/wongjj/Downloads")

df <- readRDS("analytic_sample.RDS")

# ------------------------------------------------------------------------------
# FACTOR CODING
# ------------------------------------------------------------------------------
df2 <- df %>%
  mutate(
    mvm2      = factor(mvm2,    levels = 0:1,
                       labels = c("No MVM", "Any MVM")),
    mvm3      = factor(mvm3,    levels = 0:2,
                       labels = c("None", "Low", "High")),
    ai        = factor(ai,      levels = 0:1,
                       labels = c("No AI", "Any AI")),
    ai3       = factor(ai3,     levels = 0:2,
                       labels = c("None", "Low", "High")),
    smoking   = factor(smoking, levels = 0:2,
                       labels = c("Non-smoker", "1-19/day", "20+/day")),
    parity    = factor(parity,  levels = 0:4,
                       labels = c("0", "1", "2", "3", "4+")),
    income    = factor(income,  levels = 1:7,
                       labels = c("<=$1,999", "$2,000-$3,999",
                                  "$4,000-$5,999", "$6,000-$7,999",
                                  "$8,000-$9,999", "$10,000-$14,999",
                                  ">=$15,000")),
    educ      = factor(educ,    levels = 0:3,
                       labels = c("<HS", "Some HS", "HS grad",
                                  "Some college+")),
    marital   = factor(marital, levels = 0:2,
                       labels = c("Married/CL", "Single", "Wid/Div/Sep")),
    race      = factor(race,    levels = c(1, 2, 4, 8),
                       labels = c("White", "Black", "Puerto Rican", "Other")),
    infant_sex  = factor(infant_sex,  levels = 0:1,
                         labels = c("Female", "Male")),
    dm          = factor(dm,          levels = 0:1, labels = c("No", "Yes")),
    chronic_htn = factor(chronic_htn, levels = 0:1, labels = c("No", "Yes")),
    chorangioma = factor(chorangioma, levels = 0:1, labels = c("No", "Yes")),
    mvm_villous = factor(mvm_villous, levels = 0:1, labels = c("No", "Yes")),
    mvm_vascular= factor(mvm_vascular,levels = 0:1, labels = c("No", "Yes")),
    ih1         = factor(ih1,         levels = 0:1, labels = c("No IH", "IH")),
    ih2         = factor(ih2,         levels = 0:1, labels = c("No IH", "IH"))
  )

# ------------------------------------------------------------------------------
# VARIABLE LABELS
# ------------------------------------------------------------------------------
label(df2$age)         <- "Maternal age (years)"
label(df2$bmi)         <- "Pre-pregnancy BMI (kg/m\u00b2)"
label(df2$race)        <- "Race/ethnicity"
label(df2$educ)        <- "Education"
label(df2$income)      <- "Annual income"
label(df2$marital)     <- "Marital status"
label(df2$smoking)     <- "Smoking"
label(df2$parity)      <- "Parity"
label(df2$dm)          <- "Pre-pregnancy diabetes"
label(df2$chronic_htn) <- "Chronic hypertension"
label(df2$gest_age)    <- "Gestational age (weeks)"
label(df2$birthweight) <- "Birthweight (g)"
label(df2$preterm)     <- "Preterm birth (<37 weeks)"
label(df2$SGA)         <- "Small for gestational age"
label(df2$infant_sex)  <- "Infant sex"
label(df2$chorangioma) <- "Chorangioma"
label(df2$mvm2)        <- "Maternal vascular malperfusion"
label(df2$ai)          <- "Acute inflammation"
label(df2$ih1)         <- "Infantile haemangioma at 1y (V1)"

# ------------------------------------------------------------------------------
# TABLE 1 VARIABLE LIST
# ------------------------------------------------------------------------------
tab1_var <- c("age", "bmi",
              "race", "educ", "income", "marital",
              "smoking", "parity",
              "dm", "chronic_htn",
              "gest_age", "birthweight", "preterm", "SGA",
              "infant_sex", "chorangioma",
              "ih1")

# ------------------------------------------------------------------------------
# TABLE 1 — by MVM any vs none (table1 package, pretty HTML)
# ------------------------------------------------------------------------------
aTab1_mvm <- table1(
  ~age + bmi + race + educ + income + marital +
    smoking + parity + dm + chronic_htn +
    gest_age + birthweight + preterm + SGA +
    infant_sex + chorangioma + ih1 | mvm2,
  data    = df2 %>% filter(!is.na(mvm2)),
  overall = c(right = "Overall"),
  render.continuous = c(. = "Mean (SD)", . = "Median [Q1, Q3]"))

aTab1_ai <- table1(
  ~age + bmi + race + educ + income + marital +
    smoking + parity + dm + chronic_htn +
    gest_age + birthweight + preterm + SGA +
    infant_sex + chorangioma + ih1 | ai,
  data    = df2 %>% filter(!is.na(ai)),
  overall = c(right = "Overall"),
  render.continuous = c(. = "Mean (SD)", . = "Median [Q1, Q3]"))

aTab1_ih <- table1(
  ~age + bmi + race + educ + income + marital +
    smoking + parity + dm + chronic_htn +
    gest_age + birthweight + preterm + SGA +
    infant_sex + chorangioma | ih1,
  data    = df2 %>% filter(!is.na(ih1)),
  overall = c(right = "Overall"),
  render.continuous = c(. = "Mean (SD)", . = "Median [Q1, Q3]"))

# ------------------------------------------------------------------------------
# TABLE 1 — by MVM/AI/IH (tableone package, with SMDs)
# Caniglia preference: SMD over p-values for descriptive tables
# ------------------------------------------------------------------------------
bTab1_mvm <- CreateTableOne(
  vars       = tab1_var,
  strata     = "mvm2",
  data       = df2 %>% filter(!is.na(mvm2)),
  test       = FALSE,
  includeNA  = FALSE,
  addOverall = TRUE)

bTab1_mvm_x <- print(bTab1_mvm,
                     showAllLevels = TRUE, smd = TRUE,
                     nonnormal     = "gest_age",
                     missing       = TRUE, noSpaces = TRUE,
                     printToggle   = FALSE, varLabels = TRUE)

bTab1_ai <- CreateTableOne(
  vars       = tab1_var,
  strata     = "ai",
  data       = df2 %>% filter(!is.na(ai)),
  test       = FALSE,
  includeNA  = FALSE,
  addOverall = TRUE)

bTab1_ai_x <- print(bTab1_ai,
                    showAllLevels = TRUE, smd = TRUE,
                    nonnormal     = "gest_age",
                    missing       = TRUE, noSpaces = TRUE,
                    printToggle   = FALSE, varLabels = TRUE)

bTab1_ih <- CreateTableOne(
  vars       = setdiff(tab1_var, "ih1"),
  strata     = "ih1",
  data       = df2 %>% filter(!is.na(ih1)),
  test       = FALSE,
  includeNA  = FALSE,
  addOverall = TRUE)

bTab1_ih_x <- print(bTab1_ih,
                    showAllLevels = TRUE, smd = TRUE,
                    nonnormal     = "gest_age",
                    missing       = TRUE, noSpaces = TRUE,
                    printToggle   = FALSE, varLabels = TRUE)

write.csv(bTab1_mvm_x, "table1_by_mvm_binary.csv")
write.csv(bTab1_ai_x,  "table1_by_ai_binary.csv")
write.csv(bTab1_ih_x,  "supp_table_by_ih.csv")

# ------------------------------------------------------------------------------
# CROSS-TABS (results section)
# ------------------------------------------------------------------------------
build_xtab <- function(data, exposure, exposure_label) {
  data %>%
    filter(!is.na(.data[[exposure]]), !is.na(ih1)) %>%
    group_by(.data[[exposure]]) %>%
    summarise(
      N      = n(),
      IH_n   = sum(ih1 == "IH"),
      IH_pct = round(100 * mean(ih1 == "IH"), 2),
      .groups = "drop"
    ) %>%
    rename(!!exposure_label := !!sym(exposure))
}

cat("\n========================================\n")
cat("CROSS-TABS\n")
cat("========================================\n")

# Overall
xtab_overall <- df2 %>%
  filter(!is.na(ih1)) %>%
  summarise(N = n(), IH_n = sum(ih1 == "IH"),
            IH_pct = round(100 * mean(ih1 == "IH"), 2))
cat("\nOverall IH:\n"); print(xtab_overall)

# By exposure — binary primary
xtab_mvm2  <- build_xtab(df2, "mvm2",        "MVM_any")
xtab_ai    <- build_xtab(df2, "ai",           "AI_any")

# By exposure — supplemental
xtab_mvm3  <- build_xtab(df2, "mvm3",        "MVM_grade")
xtab_ai3   <- build_xtab(df2, "ai3",         "AI_stage")
xtab_vil   <- build_xtab(df2, "mvm_villous", "MVM_villous")
xtab_ves   <- build_xtab(df2, "mvm_vascular","MVM_vascular")
xtab_chor  <- build_xtab(df2, "chorangioma", "Chorangioma")

cat("\nIH by MVM any vs none:\n");     print(xtab_mvm2)
cat("\nIH by AI any vs none:\n");      print(xtab_ai)
cat("\nIH by MVM grade (3-level):\n"); print(xtab_mvm3)
cat("\nIH by AI stage (3-level):\n");  print(xtab_ai3)
cat("\nIH by MVM villous lesions:\n"); print(xtab_vil)
cat("\nIH by MVM vascular lesions:\n");print(xtab_ves)
cat("\nIH by chorangioma:\n");         print(xtab_chor)

# MVM continuous
xtab_mvm_cont <- df2 %>%
  filter(!is.na(mvm), !is.na(ih1)) %>%
  group_by(MVM_score = mvm) %>%
  summarise(N = n(), IH_n = sum(ih1 == "IH"),
            IH_pct = round(100 * mean(ih1 == "IH"), 2), .groups = "drop")
cat("\nIH by MVM score (continuous):\n"); print(xtab_mvm_cont)

# Joint MVM x AI
xtab_joint <- df2 %>%
  filter(!is.na(mvm2), !is.na(ai), !is.na(ih1)) %>%
  mutate(joint = factor(case_when(
    mvm2 == "No MVM"  & ai == "No AI"  ~ "Neither",
    mvm2 == "No MVM"  & ai == "Any AI" ~ "AI only",
    mvm2 == "Any MVM" & ai == "No AI"  ~ "MVM only",
    mvm2 == "Any MVM" & ai == "Any AI" ~ "Both MVM + AI"
  ), levels = c("Neither","AI only","MVM only","Both MVM + AI"))) %>%
  group_by(joint) %>%
  summarise(N = n(), IH_n = sum(ih1 == "IH"),
            IH_pct = round(100 * mean(ih1 == "IH"), 2), .groups = "drop")
cat("\nIH by joint MVM x AI:\n"); print(xtab_joint)

# Save cross-tabs
write.csv(xtab_overall,  "xtab_overall.csv",       row.names = FALSE)
write.csv(xtab_mvm2,     "xtab_mvm_binary.csv",    row.names = FALSE)
write.csv(xtab_ai,       "xtab_ai_binary.csv",     row.names = FALSE)
write.csv(xtab_mvm3,     "xtab_mvm_grade.csv",     row.names = FALSE)
write.csv(xtab_ai3,      "xtab_ai_stage.csv",      row.names = FALSE)
write.csv(xtab_vil,      "xtab_mvm_villous.csv",   row.names = FALSE)
write.csv(xtab_ves,      "xtab_mvm_vascular.csv",  row.names = FALSE)
write.csv(xtab_chor,     "xtab_chorangioma.csv",   row.names = FALSE)
write.csv(xtab_mvm_cont, "xtab_mvm_continuous.csv",row.names = FALSE)
write.csv(xtab_joint,    "xtab_joint_mvm_ai.csv",  row.names = FALSE)

# ------------------------------------------------------------------------------
# DIFFERENTIAL EXCLUSION / LTFU CHECK
# MAR justification: missingness should be related to exposure (SMD > 0.1)
# ------------------------------------------------------------------------------
cat("\n========================================\n")
cat("DIFFERENTIAL EXCLUSION/LTFU BY EXPOSURE\n")
cat("========================================\n")

check_selection <- function(data, exposure) {
  data %>%
    group_by(.data[[exposure]]) %>%
    summarise(
      N                     = n(),
      n_nonsurvivor         = sum(exclusion == 1),
      pct_nonsurvivor       = round(100 * mean(exclusion == 1), 1),
      n_ltfu                = sum(exclusion == 0 & is.na(ih1)),
      pct_ltfu_of_survivors = round(100 * sum(exclusion == 0 & is.na(ih1)) /
                                      sum(exclusion == 0), 1),
      .groups = "drop"
    )
}

cat("\nNon-survivor + LTFU by MVM any vs none:\n")
print(check_selection(df2, "mvm2"))

cat("\nNon-survivor + LTFU by AI any vs none:\n")
print(check_selection(df2, "ai"))

cat("\nDone. Files saved to Downloads.\n")