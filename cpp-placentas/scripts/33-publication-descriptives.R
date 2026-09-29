# ==============================================================================
# 33-publication-descriptives.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Descriptive tables for the manuscript, rebuilt from the corrected sample
# (07). Replaces the tables made by 02 (which used 45,296 pregnancies and
# had a stale parity row, flagged by Alexa Freedman).
#
#   1. Participant flow (Figure 1): counts at each exclusion step
#   2. Table 1: characteristics overall and by any MVM and any AI
#      (45,268 pregnancies, observed data, before imputation)
#   3. Table S1: characteristics by definite IH status among the 40,700
#      infants with an observed outcome (Overall column = 40,700, per
#      Alexa Freedman's comment)
#   4. Table S2: study site description (pregnancies, follow-up, MVM, AI,
#      IH) in the corrected sample
#   5. Missingness of covariates (for the Methods)
#
# Format: percentages before counts, e.g. 29.0% (13,119). Continuous
# variables as mean (SD). Missing values shown as their own row.
#
# Outputs (outputs/): flow_counts.csv, table1_corrected.csv,
#   tableS1_by_ih.csv, tableS2_site.csv, missingness.csv
# ==============================================================================

library(tidyverse)
library(haven)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")

clean_numeric <- function(x) {
  x <- as.character(x); x[x %in% c(".u", ".U", "u", "U")] <- NA; as.numeric(x)
}
n_ <- function(x) format(x, big.mark = ",", trim = TRUE)
pct_n <- function(k, n) sprintf("%.1f%% (%s)", 100 * k / n, n_(k))

# ------------------------------------------------------------------------------
# 1. PARTICIPANT FLOW (same steps as 07)
# ------------------------------------------------------------------------------
ncpp <- read_sas(file.path(in_dir, "ncppbasa.sas7bdat"),
                 col_select = c(MOMID, PREGID, CHILDID, C10))
p1 <- read_sas(file.path(in_dir, "path1all.sas7bdat"),
               col_select = c(MOMID, PREGID, CHILDID))
p2 <- read_sas(file.path(in_dir, "path2all.sas7bdat"),
               col_select = c(MOMID, PREGID, CHILDID))
mvm <- read_csv(file.path(in_dir, "CPP_dataset_021226.csv"), show_col_types = FALSE) %>%
  rename(PREGID = NINDB_NUMBER) %>% group_by(PREGID) %>%
  slice_max(MVM_SCORE, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(PREGID, MVM_SCORE, AI_DI)

step0 <- ncpp
path  <- p1 %>% inner_join(p2, by = c("MOMID", "PREGID", "CHILDID"))
step1 <- step0 %>% inner_join(path, by = c("MOMID", "PREGID", "CHILDID"))
dup   <- step1 %>% count(PREGID) %>% filter(n > 1) %>% pull(PREGID)
step2 <- step1 %>% filter(!PREGID %in% dup, clean_numeric(C10) == 1)
step3 <- step2 %>% left_join(mvm, by = "PREGID") %>% filter(!is.na(MVM_SCORE))
step4 <- step3 %>% filter(!is.na(AI_DI))

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS"))
stopifnot(nrow(step4) == nrow(d))

flow <- tibble(
  step = c("Pregnancies (infants) in the CPP",
           "Excluded: no placental pathology record",
           "Infants with a placental pathology record",
           "Excluded: multiple gestation",
           "Excluded: missing MVM score",
           "Excluded: missing AI",
           "Analytic sample (Table 1)",
           "Excluded: died before the one-year examination",
           "Excluded: alive but no one-year IH examination",
           "Infants with an observed IH outcome (Table 2)",
           "Definite IH cases"),
  n = c(nrow(step0), nrow(step0) - nrow(step1), nrow(step1),
        nrow(step1) - nrow(step2), nrow(step2) - nrow(step3),
        nrow(step3) - nrow(step4), nrow(d),
        sum(d$exclusion == 1), sum(d$exclusion == 0 & is.na(d$ih1)),
        sum(!is.na(d$ih1) & d$exclusion == 0), sum(d$ih1 == 1, na.rm = TRUE)))
rm(ncpp, p1, p2, path, dup, step0, step1, step2, step3, step4); invisible(gc())

# ------------------------------------------------------------------------------
# TABLE HELPERS
# ------------------------------------------------------------------------------
d <- d %>% mutate(
  age_cat   = cut(age, c(-Inf, 19, 24, 29, 34, Inf),
                  labels = c("<20", "20-24", "25-29", "30-34", ">=35")),
  race      = factor(race, c(1, 2, 4, 8), c("White", "Black", "Puerto Rican", "Other")),
  educ      = factor(educ, 0:3, c("Less than high school", "Some high school",
                                  "High school graduate", "Some college or more")),
  income    = factor(income, 1:7, c("$1,999 or less", "$2,000 to $3,999",
                                    "$4,000 to $5,999", "$6,000 to $7,999",
                                    "$8,000 to $9,999", "$10,000 to $14,999",
                                    "$15,000 or more")),
  marital   = factor(marital, 0:2, c("Married or common law", "Single",
                                     "Widowed, divorced, or separated")),
  smoking   = factor(smoking, 0:2, c("Non-smoker", "1 to 19 cigarettes per day",
                                     "20 or more cigarettes per day")),
  parity    = factor(parity, 0:4, c("0", "1", "2", "3", "4 or more")),
  dm        = factor(dm, 0:1, c("No", "Yes")),
  chronic_htn = factor(chronic_htn, 0:1, c("No", "Yes")),
  infant_sex = factor(infant_sex, 0:1, c("Female", "Male")),
  preterm   = factor(preterm, 0:1, c("No", "Yes")),
  SGA       = factor(SGA, 0:1, c("No", "Yes")),
  mvm_any   = factor(as.integer(mvm2 == 1), 0:1, c("None", "Any")),
  ai_any    = factor(as.integer(ai == 1), 0:1, c("None", "Any")),
  mvm_grade = factor(mvm3, 0:2, c("None", "Low (score 1 to 2)", "High (score 3 or more)")),
  ai_stage  = factor(ai3, 0:2, c("None", "Low", "High")),
  joint     = factor(case_when(mvm2 == 0 & ai == 0 ~ "Neither",
                               mvm2 == 1 & ai == 0 ~ "MVM only",
                               mvm2 == 0 & ai == 1 ~ "AI only",
                               TRUE ~ "Both MVM and AI"),
                     c("Neither", "MVM only", "AI only", "Both MVM and AI")))

cont_row <- function(dd, v, digits = 1) {
  x <- dd[[v]]
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}
cat_rows <- function(dd, v) {
  x <- dd[[v]]; n_obs <- sum(!is.na(x))
  out <- sapply(levels(x), function(l) pct_n(sum(x == l, na.rm = TRUE), n_obs))
  if (any(is.na(x))) out <- c(out, Missing = n_(sum(is.na(x))))
  out
}

vars <- list(
  list("Maternal age, years, mean (SD)", "age", "cont"),
  list("Pre-pregnancy BMI, kg/m², mean (SD)", "bmi", "cont"),
  list("Race and ethnicity", "race", "cat"),
  list("Education", "educ", "cat"),
  list("Family income, per year", "income", "cat"),
  list("Marital status", "marital", "cat"),
  list("Smoking during pregnancy", "smoking", "cat"),
  list("Parity", "parity", "cat"),
  list("Pre-pregnancy diabetes", "dm", "cat"),
  list("Chronic hypertension", "chronic_htn", "cat"),
  list("Infant sex", "infant_sex", "cat"),
  list("Perinatal characteristics (not adjusted for)", NA, "header"),
  list("Gestational age at delivery, weeks, mean (SD)", "gest_age", "cont"),
  list("Preterm birth (<37 weeks)", "preterm", "cat"),
  list("Birthweight, g, mean (SD)", "birthweight", "cont0"),
  list("Small for gestational age", "SGA", "cat"))

build_table <- function(groups) {
  map_dfr(vars, function(v) {
    lab <- v[[1]]; var <- v[[2]]; type <- v[[3]]
    if (type == "header") return(tibble(Characteristic = lab))
    if (type %in% c("cont", "cont0")) {
      vals <- map_chr(groups, ~ cont_row(.x, var, ifelse(type == "cont0", 0, 1)))
      miss <- map_int(groups, ~ sum(is.na(.x[[var]])))
      out <- tibble(Characteristic = lab, !!!setNames(as.list(vals), names(groups)))
      if (any(miss > 0))
        out <- bind_rows(out, tibble(Characteristic = "   Missing",
                                     !!!setNames(as.list(n_(miss)), names(groups))))
      return(out)
    }
    rows <- map(groups, ~ cat_rows(.x, var))
    lv <- unique(unlist(map(rows, names)))
    bind_rows(tibble(Characteristic = lab),
              tibble(Characteristic = paste0("   ", lv),
                     !!!map(rows, ~ unname(ifelse(lv %in% names(.x), .x[lv],
                                                  ifelse(lv == "Missing", "0", ""))))))
  }) %>% mutate(across(everything(), ~ replace_na(.x, "")))
}

# ------------------------------------------------------------------------------
# 2. TABLE 1 (45,268 pregnancies)
# ------------------------------------------------------------------------------
g1 <- list(Overall = d,
           `No MVM` = filter(d, mvm_any == "None"), `Any MVM` = filter(d, mvm_any == "Any"),
           `No AI` = filter(d, ai_any == "None"), `Any AI` = filter(d, ai_any == "Any"))
t1 <- build_table(g1)
hdr1 <- map_chr(g1, ~ sprintf("n = %s", n_(nrow(.x))))
t1 <- bind_rows(tibble(Characteristic = "", !!!as.list(hdr1)), t1)

# ------------------------------------------------------------------------------
# 3. TABLE S1 (40,700 infants with observed outcome, by definite IH)
# ------------------------------------------------------------------------------
obs <- d %>% filter(!is.na(ih1), exclusion == 0)
stopifnot(nrow(obs) == 40700, sum(obs$ih1) == 700)
g2 <- list(Overall = obs, `No IH` = filter(obs, ih1 == 0), `Definite IH` = filter(obs, ih1 == 1))
vars_s1 <- vars
vars <- c(list(list("Placental pathology", NA, "header"),
               list("Maternal vascular malperfusion", "mvm_any", "cat"),
               list("MVM grade", "mvm_grade", "cat"),
               list("MVM score, mean (SD)", "mvm", "cont"),
               list("Acute inflammation", "ai_any", "cat"),
               list("AI stage", "ai_stage", "cat"),
               list("Joint MVM and AI", "joint", "cat"),
               list("Maternal and infant characteristics", NA, "header")), vars_s1)
s1 <- build_table(g2)
hdr2 <- map_chr(g2, ~ sprintf("n = %s", n_(nrow(.x))))
s1 <- bind_rows(tibble(Characteristic = "", !!!as.list(hdr2)), s1)
vars <- vars_s1

# ------------------------------------------------------------------------------
# 4. TABLE S2: STUDY SITE (corrected sample)
# ------------------------------------------------------------------------------
site_names <- c("5" = "Boston", "10" = "Buffalo", "15" = "New Orleans",
                "31" = "New York (Columbia)", "37" = "Baltimore", "45" = "Virginia",
                "50" = "Minnesota", "55" = "New York Medical College", "60" = "Oregon",
                "66" = "Pennsylvania", "71" = "Providence", "82" = "Tennessee")
site_row <- function(dd, lab) {
  o <- dd %>% filter(!is.na(ih1), exclusion == 0)
  tibble(Site = lab,
         Pregnancies = n_(nrow(dd)),
         `Any MVM` = pct_n(sum(dd$mvm2 == 1), nrow(dd)),
         `Any AI` = pct_n(sum(dd$ai == 1), nrow(dd)),
         `Observed at one year` = pct_n(nrow(o), nrow(dd)),
         `Definite IH, % (cases/N)` = sprintf("%.1f%% (%s/%s)", 100 * mean(o$ih1),
                                              n_(sum(o$ih1)), n_(nrow(o))))
}
s2 <- bind_rows(
  map_dfr(names(site_names), ~ site_row(filter(d, as.character(site) == .x),
                                        paste(.x, site_names[.x]))),
  site_row(d, "All sites"))

# ------------------------------------------------------------------------------
# 5. MISSINGNESS
# ------------------------------------------------------------------------------
miss <- tibble(variable = c("bmi", "income", "educ", "parity", "smoking",
                            "dm", "chronic_htn", "marital")) %>%
  mutate(n_missing = map_int(variable, ~ sum(is.na(d[[.x]]))),
         pct = sprintf("%.1f%%", 100 * n_missing / nrow(d)))

print(flow, n = Inf); print(t1, n = Inf, width = Inf); print(s2, width = Inf); print(miss)
write_csv(flow, file.path(out_dir, "flow_counts.csv"))
write_csv(t1, file.path(out_dir, "table1_corrected.csv"))
write_csv(s1, file.path(out_dir, "tableS1_by_ih.csv"))
write_csv(s2, file.path(out_dir, "tableS2_site.csv"))
write_csv(miss, file.path(out_dir, "missingness.csv"))
cat("\nSaved: flow_counts, table1_corrected, tableS1_by_ih, tableS2_site, missingness\n")
