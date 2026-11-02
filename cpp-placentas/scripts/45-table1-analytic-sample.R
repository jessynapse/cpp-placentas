# ==============================================================================
# 45-table1-analytic-sample.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Table 1 on the ANALYTIC SAMPLE: the 40,700 singleton infants examined for
# IH at 12 months (the sample used in every outcome model), overall and by
# any MVM and any AI. Script 33 built Table 1 on the 45,268 eligible
# pregnancies. The PPE abstract and Results now describe the 40,700, so
# Table 1 matches them (as in Freedman et al. 2025, who described the
# sample with the outcome measured).
#
# Data build: reproduces 07 (corrected pipeline), including the inner join
# of path1all and path2all that defines a placental pathology record.
# path2all alone gives 47,733 eligible pregnancies, not 45,268, so path1all
# is required. The script stops unless the eligible cohort is exactly
# 45,268, the analytic sample 40,700 and definite IH 700.
#
# Format: percentages before counts, e.g. 27.8% (11,324). Continuous
# variables as mean (SD). Missing values shown as their own row.
#
# Inputs: data/ (or /Users/wongjj/Downloads): ncppbasa.sas7bdat,
#   path1all.sas7bdat, path2all.sas7bdat, cpp_diabetes_20231026.sas7bdat,
#   cpp_htn_20231031.sas7bdat, CPP_dataset_021226.csv
# Output: outputs/table1_analytic.csv
# ==============================================================================

library(tidyverse)
library(haven)

in_dir <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
if (!file.exists(file.path(in_dir, "ncppbasa.sas7bdat")) && dir.exists("../data")) in_dir <- "../data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"

clean_numeric <- function(x) {
  x <- as.character(x); x[x %in% c(".u", ".U", "u", "U")] <- NA; as.numeric(x)
}
n_ <- function(x) format(x, big.mark = ",", trim = TRUE)
pct_n <- function(k, n) sprintf("%.1f%% (%s)", 100 * k / n, n_(k))

# ------------------------------------------------------------------------------
# 1. DATA BUILD (as in 07)
# ------------------------------------------------------------------------------
path1 <- read_sas(file.path(in_dir, "path1all.sas7bdat"),
                  col_select = c(MOMID, PREGID, CHILDID))
path2 <- read_sas(file.path(in_dir, "path2all.sas7bdat"),
                  col_select = c(MOMID, PREGID, CHILDID))
ncpp <- read_sas(file.path(in_dir, "ncppbasa.sas7bdat"),
                 col_select = c(MOMID, PREGID, CHILDID, C983, C1095, C1101,
                                C31, C83, C81, C52, C50, C10, C275, C272,
                                C36, C554, C303, C1092, SITE))
dm  <- read_sas(file.path(in_dir, "cpp_diabetes_20231026.sas7bdat")) %>%
  select(MOMID, PREGID, DM)
htn <- read_sas(file.path(in_dir, "cpp_htn_20231031.sas7bdat")) %>%
  select(MOMID, PREGID, HTN_PREG) %>%
  mutate(chronic_htn = if_else(is.na(HTN_PREG), NA_integer_,
                               as.integer(HTN_PREG %in% c(1, 4))))
mvm <- read_csv(file.path(in_dir, "CPP_dataset_021226.csv"), show_col_types = FALSE) %>%
  rename(PREGID = NINDB_NUMBER) %>% group_by(PREGID) %>%
  slice_max(MVM_SCORE, n = 1, with_ties = FALSE) %>% ungroup()

d <- path1 %>%
  inner_join(path2, by = c("MOMID", "PREGID", "CHILDID")) %>%
  inner_join(ncpp, by = c("MOMID", "PREGID", "CHILDID")) %>%
  left_join(dm,  by = c("MOMID", "PREGID")) %>%
  left_join(htn, by = c("MOMID", "PREGID")) %>%
  left_join(mvm, by = "PREGID") %>%
  rename(birthweight = C1095, gest_age = C1101) %>%
  mutate(IH_raw    = clean_numeric(C983),
         exclusion = ifelse(C1092 == 0 | (C1092 >= 51 & C1092 <= 58), 0, 1),
         ih1       = case_when(exclusion == 1 ~ NA_real_,
                               IH_raw == 2    ~ 0,
                               TRUE           ~ IH_raw))

twins <- d %>% count(PREGID) %>% filter(n > 1) %>% pull(PREGID)
d <- d %>% filter(!PREGID %in% twins, !is.na(MVM_SCORE),
                  clean_numeric(C10) == 1, !is.na(AI_DI))

d <- d %>% mutate(
  age     = clean_numeric(C31),
  bmi     = ifelse(clean_numeric(C81) * 0.0254 > 0,
                   (clean_numeric(C83) * 0.453592) / (clean_numeric(C81) * 0.0254)^2, NA),
  smoking = case_when(clean_numeric(C52) == 0 ~ 0,
                      clean_numeric(C52) >= 1 & clean_numeric(C52) <= 19 ~ 1,
                      clean_numeric(C52) >= 20 ~ 2),
  parity  = pmin(case_when(na_tag(C50) == "p" ~ 0,
                           na_tag(C50) == "u" ~ NA_real_,
                           is.na(C50)         ~ NA_real_,
                           TRUE               ~ as.numeric(C50)), 4),
  income  = clean_numeric(C275),
  educ    = case_when(clean_numeric(C272) %in% 0:2 ~ 0, clean_numeric(C272) == 3 ~ 1,
                      clean_numeric(C272) == 4 ~ 2, clean_numeric(C272) %in% 5:6 ~ 3),
  marital = case_when(clean_numeric(C36) %in% c(2, 3) ~ 0, clean_numeric(C36) == 1 ~ 1,
                      clean_numeric(C36) %in% c(4, 5, 6) ~ 2),
  infant_sex = case_when(clean_numeric(C554) == 1 ~ 1, clean_numeric(C554) == 2 ~ 0),
  race    = case_when(clean_numeric(C303) == 3 ~ 8, TRUE ~ clean_numeric(C303)),
  mvm2    = as.integer(MVM_SCORE_CAT > 0),
  ai      = AI_DI,
  preterm = as.integer(gest_age < 37)) %>%
  group_by(gest_age, infant_sex) %>%
  mutate(SGA = as.integer((birthweight - mean(birthweight, na.rm = TRUE)) /
                            sd(birthweight, na.rm = TRUE) < -1.28)) %>%
  ungroup()

obs <- d %>% filter(!is.na(ih1), exclusion == 0)
cat("Eligible cohort:", nrow(d), " Analytic sample:", nrow(obs),
    " Definite IH:", sum(obs$ih1), " Any MVM:", sum(obs$mvm2), " Any AI:", sum(obs$ai), "\n")
stopifnot(nrow(d) == 45268, nrow(obs) == 40700, sum(obs$ih1) == 700,
          sum(obs$mvm2) == 11324, sum(obs$ai) == 8880)

# ------------------------------------------------------------------------------
# 2. TABLE 1 (40,700 infants)
# ------------------------------------------------------------------------------
obs <- obs %>% mutate(
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
  dm        = factor(DM, 0:1, c("No", "Yes")),
  chronic_htn = factor(chronic_htn, 0:1, c("No", "Yes")),
  infant_sex = factor(infant_sex, 0:1, c("Female", "Male")),
  preterm   = factor(preterm, 0:1, c("No", "Yes")),
  SGA       = factor(SGA, 0:1, c("No", "Yes")),
  mvm_any   = factor(mvm2, 0:1, c("None", "Any")),
  ai_any    = factor(ai, 0:1, c("None", "Any")))

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

g1 <- list(Overall = obs,
           `No MVM` = filter(obs, mvm_any == "None"), `Any MVM` = filter(obs, mvm_any == "Any"),
           `No AI` = filter(obs, ai_any == "None"), `Any AI` = filter(obs, ai_any == "Any"))
t1 <- build_table(g1)
hdr <- map_chr(g1, ~ sprintf("n = %s", n_(nrow(.x))))
t1 <- bind_rows(tibble(Characteristic = "", !!!as.list(hdr)), t1)

dir.create(out_dir, showWarnings = FALSE)
write_csv(t1, file.path(out_dir, "table1_analytic.csv"))
cat("Saved:", file.path(out_dir, "table1_analytic.csv"), "\n")
