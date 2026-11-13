# ==============================================================================
# 68-table1-joint-exposure.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Table 1 by the four-category joint exposure (now the main exposure):
# neither (reference), maternal vascular malperfusion only, acute
# inflammation only, both. Same data build, birthweight cleaning (under
# 300 g or over 6,000 g set to missing), rows, labels and format as 67.
#
# SMD: one column, the LARGEST of the three SMDs comparing each exposure
# category with neither (Yang and Dalton 2012 for categorical variables,
# missing excluded). The three separate SMDs are saved alongside in
# outputs/table1_joint_smd_detail.csv.
#
# Outputs: outputs/table1_joint.csv, outputs/table1_joint_smd_detail.csv
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


# ------------------------------------------------------------------------------
# 1b. SET IMPLAUSIBLE BIRTHWEIGHT TO MISSING
# ------------------------------------------------------------------------------
d <- d %>% mutate(bw_raw = clean_numeric(birthweight), examined_ = !is.na(ih1) & exclusion == 0,
                  rule_bw = !is.na(bw_raw) & (bw_raw < 300 | bw_raw > 6000),
                  birthweight = ifelse(rule_bw, NA, bw_raw))
clean_counts <- tibble(
  rule = c("Birthweight under 300 g or over 6,000 g", "Birthweight missing after cleaning"),
  all_45268 = c(sum(d$rule_bw), sum(is.na(d$birthweight))),
  examined_40700 = c(sum(d$rule_bw & d$examined_), sum(is.na(d$birthweight) & d$examined_)))
print(clean_counts)
d <- d %>% select(-bw_raw, -examined_, -rule_bw)

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
# SGA needs gestational age: missing GA was being grouped as its own "NA" stratum in 45
d <- d %>% mutate(SGA = ifelse(is.na(gest_age) | is.na(birthweight), NA_integer_, SGA))

obs <- d %>% filter(!is.na(ih1), exclusion == 0)
cat("Eligible cohort:", nrow(d), " Analytic sample:", nrow(obs),
    " Definite IH:", sum(obs$ih1), " Any MVM:", sum(obs$mvm2), " Any AI:", sum(obs$ai), "\n")
stopifnot(nrow(d) == 45268, nrow(obs) == 40700, sum(obs$ih1) == 700,
          sum(obs$mvm2) == 11324, sum(obs$ai) == 8880)

# ------------------------------------------------------------------------------
# 2. TABLE 1 (40,700 infants), laid out like the ECHO Cycle 2 Table 1:
#    No. (%) with count first, Missing as its own row (% of column total),
#    continuous variables as median (IQR) and range plus categories,
#    binary variables on one row, and a standardised mean difference (SMD)
#    for each characteristic (exposed vs unexposed, missing excluded).
# ------------------------------------------------------------------------------
obs <- obs %>% mutate(
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
  infant_sex = factor(infant_sex, 0:1, c("Female", "Male")))

obs <- obs %>% mutate(joint = case_when(mvm2 == 0 & ai == 0 ~ "Neither", mvm2 == 1 & ai == 0 ~ "MVM only",
                                        mvm2 == 0 & ai == 1 ~ "AI only", TRUE ~ "Both"))
grp <- list(Total = obs, Neither = filter(obs, joint == "Neither"), `MVM only` = filter(obs, joint == "MVM only"),
            `AI only` = filter(obs, joint == "AI only"), Both = filter(obs, joint == "Both"))
stopifnot(map_int(grp, nrow) == c(40700, 23391, 8429, 5985, 2895))

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
smd_detail <- list()
smd_max <- function(var, type, label) {
  f <- if (type == "cont") smd_cont else smd_cat
  s <- map_dbl(c("MVM only", "AI only", "Both"), ~ f(grp[[.x]][[var]], grp$Neither[[var]]))
  smd_detail[[length(smd_detail) + 1]] <<- tibble(Characteristic = label, `MVM only vs neither` = round(s[1], 2),
                                                 `AI only vs neither` = round(s[2], 2), `Both vs neither` = round(s[3], 2))
  sprintf("%.2f", max(s))
}
mk <- function(label, vals = rep("", 5), smd = "") {
  tibble(Characteristic = label, Total = vals[1], Neither = vals[2], `MVM only` = vals[3],
         `AI only` = vals[4], Both = vals[5], SMD = smd)
}
miss_row <- function(var) {
  k <- map_int(grp, ~ sum(is.na(.x[[var]])))
  if (k[1] == 0) return(NULL)
  mk("   Missing", unname(map2_chr(k, map_int(grp, nrow), np)))
}
cont_block <- function(label, var, dg, catvar = NULL) {
  med <- map_chr(grp, ~ { q <- quantile(.x[[var]], c(.5, .25, .75), na.rm = TRUE)
                          sprintf("%s (%s-%s)", fmt(q[1], dg), fmt(q[2], dg), fmt(q[3], dg)) })
  rng <- map_chr(grp, ~ sprintf("%s-%s", fmt(min(.x[[var]], na.rm = TRUE), dg),
                                fmt(max(.x[[var]], na.rm = TRUE), dg)))
  out <- bind_rows(mk(label, smd = smd_max(var, "cont", label)),
                   mk("   Median (IQR)", unname(med)), mk("   Range", unname(rng)))
  if (!is.null(catvar)) out <- bind_rows(out, cat_levels(catvar))
  bind_rows(out, miss_row(var))
}
cat_levels <- function(var) {
  lv <- levels(obs[[var]])
  map_dfr(lv, function(l) mk(paste0("   ", l), unname(map_chr(grp, ~ {
    x <- .x[[var]]; np(sum(x == l, na.rm = TRUE), sum(!is.na(x))) }))))
}
cat_block <- function(label, var) bind_rows(mk(label, smd = smd_max(var, "cat", label)), cat_levels(var), miss_row(var))
bin_block <- function(label, var) {
  vals <- map_chr(grp, ~ { x <- .x[[var]]; np(sum(x == 1, na.rm = TRUE), sum(!is.na(x))) })
  grp <<- map(grp, ~ mutate(.x, !!paste0(var, "_f") := factor(.x[[var]], 0:1)))
  bind_rows(mk(label, unname(vals), smd_max(paste0(var, "_f"), "cat", label)), miss_row(var))
}

t1 <- bind_rows(
  mk("", unname(map_chr(grp, ~ sprintf("n = %s", n_(nrow(.x)))))),
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
t1 <- t1 %>% filter(!(Characteristic == "   Missing" & lag(Characteristic) == "   Missing"))

dir.create(out_dir, showWarnings = FALSE)
write_csv(t1, file.path(out_dir, "table1_joint.csv"))
write_csv(bind_rows(smd_detail), file.path(out_dir, "table1_joint_smd_detail.csv"))
cat("Saved: table1_joint.csv, table1_joint_smd_detail.csv\n")
