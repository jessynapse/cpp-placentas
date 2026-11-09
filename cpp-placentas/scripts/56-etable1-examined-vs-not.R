# ==============================================================================
# 56-etable1-examined-vs-not.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# eTable 1 (examined vs not examined at one year) rebuilt with EXACTLY the
# same variables, labels, order, summaries and format as Table 1 (45), so
# the main and supplementary characteristics tables match:
#   - same sections and variable order as Table 1
#   - No. (%) with count first, Missing rows as % of column total
#   - continuous variables as median (IQR) and range plus the same categories
#   - binary variables on one row
#   - SMD for examined vs all not examined (missing excluded)
# Adds a placental pathology section (any MVM, any AI, joint exposure),
# since these are the stratifiers in Table 1.
# Data build and row builders copied from 45 (not edited).
#
# Columns: Examined (40,700), Not examined total (4,568), Died before one
# year (1,834), Alive without one-year examination (2,734).
#
# Inputs: as 45. Output: outputs/etable1_examined_vs_not.csv
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
pct_n <- function(k, n) ifelse(k > 0 & 100 * k / n < 0.05, sprintf("<0.1%% (%s)", n_(k)),
                               sprintf("%.1f%% (%s)", 100 * k / n, n_(k)))

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
  parity_raw = case_when(na_tag(C50) == "p" ~ 0,
                         na_tag(C50) == "u" ~ NA_real_,
                         is.na(C50)         ~ NA_real_,
                         TRUE               ~ as.numeric(C50)),
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
# 2. eTABLE 1 (45,268 pregnancies), same layout as Table 1
# ------------------------------------------------------------------------------
d <- d %>% mutate(
  age_cat = cut(age, c(-Inf, 19.999, 24.999, 29.999, 34.999, Inf),
                labels = c("<20", "20-24", "25-29", "30-34", ">=35")),
  race    = factor(race, c(1, 2, 4, 8), c("White", "Black", "Puerto Rican", "Other")),
  educ    = factor(educ, 0:3, c("Less Than High School", "Some High School",
                                "High School Graduate", "Some College or More")),
  income  = factor(income, 1:7, c("<=$1,999", "$2,000-$3,999", "$4,000-$5,999",
                                  "$6,000-$7,999", "$8,000-$9,999", "$10,000-$14,999",
                                  ">=$15,000")),
  marital = factor(marital, 0:2, c("Married or Common Law", "Single",
                                   "Widowed, Divorced, or Separated")),
  parity_cat = factor(parity, 0:4, c("0", "1", "2", "3", ">=4")),
  bmi_cat = cut(bmi, c(-Inf, 18.5, 25, 30, Inf), right = FALSE,
                labels = c("Underweight (<18.5)", "Normal Weight (18.5-24.9)",
                           "Overweight (25.0-29.9)", "Obesity (>=30.0)")),
  smoking = factor(smoking, 0:2, c("Non-Smoker", "1-19 Cigarettes/d", ">=20 Cigarettes/d")),
  infant_sex = factor(infant_sex, 0:1, c("Female", "Male")),
  mvm_f   = factor(mvm2, 1:0, c("Any", "None")),
  ai_f    = factor(ai, 1:0, c("Any", "None")),
  joint   = factor(case_when(mvm2 == 0 & ai == 0 ~ "Neither", mvm2 == 1 & ai == 0 ~ "Maternal Vascular Malperfusion Only",
                             mvm2 == 0 & ai == 1 ~ "Acute Inflammation Only", TRUE ~ "Both"),
                   levels = c("Neither", "Maternal Vascular Malperfusion Only", "Acute Inflammation Only", "Both")),
  examined = !is.na(ih1) & exclusion == 0)

grp <- list(Examined = filter(d, examined), `Not examined` = filter(d, !examined),
            Died = filter(d, exclusion == 1), `Alive, no examination` = filter(d, exclusion == 0, is.na(ih1)))
stopifnot(map_int(grp, nrow) == c(40700, 4568, 1834, 2734))

np <- function(k, n) {
  p <- 100 * k / n
  sprintf("%s (%s)", n_(k), ifelse(k > 0 & p < 0.05, "<0.1", sprintf("%.1f", p)))
}
fmt <- function(x, dg) formatC(x, format = "f", digits = dg, big.mark = ",")
smd_cont <- function(x1, x2) {
  abs(mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)) /
    sqrt((var(x1, na.rm = TRUE) + var(x2, na.rm = TRUE)) / 2)
}
smd_cat <- function(x1, x2) {           # Yang and Dalton (2012), missing excluded
  lv <- levels(x1)
  p1 <- prop.table(table(factor(x1[!is.na(x1)], lv)))
  p2 <- prop.table(table(factor(x2[!is.na(x2)], lv)))
  if (length(lv) == 2) {
    a <- p1[2]; b <- p2[2]
    return(as.numeric(abs(a - b) / sqrt((a * (1 - a) + b * (1 - b)) / 2)))
  }
  dd <- (p1 - p2)[-1]; q1 <- p1[-1]; q2 <- p2[-1]
  S <- (-(outer(q1, q1) + outer(q2, q2)) / 2)
  diag(S) <- (q1 * (1 - q1) + q2 * (1 - q2)) / 2
  as.numeric(sqrt(t(dd) %*% solve(S) %*% dd))
}
# one SMD: examined vs all not examined
smd_one <- function(var, type) {
  f <- if (type == "cont") smd_cont else smd_cat
  sprintf("%.2f", f(grp$Examined[[var]], grp$`Not examined`[[var]]))
}
mk <- function(label, vals = rep("", 4), smd = "") {
  tibble(Characteristic = label, Examined = vals[1], `Not examined, total` = vals[2],
         `Died before one year` = vals[3], `Alive, no one-year examination` = vals[4], SMD = smd)
}
miss_row <- function(var) {
  k <- map_int(grp, ~ sum(is.na(.x[[var]])))
  if (all(k == 0)) return(NULL)
  mk("   Missing", unname(map2_chr(k, map_int(grp, nrow), np)))
}
cont_block <- function(label, var, dg, catvar = NULL) {
  med <- map_chr(grp, ~ { q <- quantile(.x[[var]], c(.5, .25, .75), na.rm = TRUE)
                          sprintf("%s (%s-%s)", fmt(q[1], dg), fmt(q[2], dg), fmt(q[3], dg)) })
  rng <- map_chr(grp, ~ sprintf("%s-%s", fmt(min(.x[[var]], na.rm = TRUE), dg),
                                fmt(max(.x[[var]], na.rm = TRUE), dg)))
  out <- bind_rows(mk(label, smd = smd_one(var, "cont")),
                   mk("   Median (IQR)", unname(med)), mk("   Range", unname(rng)))
  if (!is.null(catvar)) out <- bind_rows(out, cat_levels(catvar))
  bind_rows(out, miss_row(var))
}
cat_levels <- function(var) {
  lv <- levels(d[[var]])
  map_dfr(lv, function(l) mk(paste0("   ", l), unname(map_chr(grp, ~ {
    x <- .x[[var]]; np(sum(x == l, na.rm = TRUE), sum(!is.na(x))) }))))
}
cat_block <- function(label, var, smd_var = var) {
  bind_rows(mk(label, smd = smd_one(smd_var, "cat")), cat_levels(var), miss_row(var))
}
bin_block <- function(label, var) {
  vals <- map_chr(grp, ~ { x <- .x[[var]]; np(sum(x == 1, na.rm = TRUE), sum(!is.na(x))) })
  grp <<- map(grp, ~ mutate(.x, !!paste0(var, "_f") := factor(.x[[var]], 0:1)))
  bind_rows(mk(label, unname(vals), smd_one(paste0(var, "_f"), "cat")), miss_row(var))
}

et1 <- bind_rows(
  mk("", unname(map_chr(grp, ~ sprintf("n = %s", n_(nrow(.x)))))),
  mk("## Placental Pathology"),
  cat_block("Maternal Vascular Malperfusion", "mvm_f"),
  cat_block("Acute Inflammation", "ai_f"),
  cat_block("Joint Exposure", "joint"),
  mk("## Maternal Sociodemographic Characteristics"),
  cont_block("Age at Delivery, y", "age", 0, "age_cat"),
  cat_block("Race and Ethnicity", "race"),
  cat_block("Educational Attainment", "educ"),
  cat_block("Annual Family Income", "income"),
  cat_block("Marital Status", "marital"),
  mk("## Maternal Health and Pregnancy Exposures"),
  cont_block("Number of Prior Births", "parity_raw", 0, "parity_cat"),
  cont_block("Pre-Pregnancy BMI", "bmi", 1, "bmi_cat"),
  bin_block("Pre-Pregnancy Diabetes", "DM"),
  bin_block("Chronic Hypertension", "chronic_htn"),
  cat_block("Smoking During Pregnancy", "smoking"),
  mk("## Infant Characteristics"),
  cat_block("Sex", "infant_sex"),
  cont_block("Gestational Age at Delivery, wk^a", "gest_age", 0),
  bin_block("Preterm Birth (<37 wk)^a", "preterm"),
  cont_block("Birth Weight, g^a", "birthweight", 0),
  bin_block("Small for Gestational Age^a", "SGA"))

# parity categories come from the capped variable, so their Missing row
# duplicates the raw-parity Missing row: keep one (as in 45)
et1 <- et1 %>% filter(!(Characteristic == "   Missing" & lag(Characteristic) == "   Missing"))

dir.create(out_dir, showWarnings = FALSE)
write_csv(et1, file.path(out_dir, "etable1_examined_vs_not.csv"))
print(et1, n = Inf, width = Inf)
cat("Saved:", file.path(out_dir, "etable1_examined_vs_not.csv"), "\n")
