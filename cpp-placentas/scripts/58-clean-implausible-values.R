# ==============================================================================
# 58-clean-implausible-values.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Sets implausible or miscoded values to missing in the corrected analytic
# sample (07), before the full rerun (59 imputation -> 60 IPCW -> 61 models).
# 07 itself is not edited. Covariates set to missing are then imputed.
#
# Covariates (affect the models):
#   1. Maternal age under 13 years                       -> age missing
#   2. Height under 122 cm or over 213 cm                -> BMI missing
#   3. Pre-pregnancy weight under 30 kg                  -> BMI missing
#   4. BMI under 13 or over 60 kg/m2                     -> BMI missing
#   5. Cigarettes per day coded 70 or 80                 -> smoking missing
#      These are codes, not counts: term infants of these mothers weigh
#      3,310 and 3,329 g on average, the same as non-smokers (3,303 g) and
#      above every true smoking category (3,117 to 3,176 g). 07 had placed
#      them in the 20+/day group.
# Descriptive only (not covariates or imputation predictors):
#   6. Gestational age outside 20 to 44 completed weeks  -> missing
#   7. Birthweight under 300 g or over 6,000 g           -> missing
#   8. Birthweight implausible for gestational age       -> gestational age
#      missing (GA under 28 wk with BW over 2,500 g, or GA 37 wk or more
#      with BW under 1,000 g)
#   Preterm birth, preterm3, SGA and the BW-for-GA z-score are re-derived
#   from the cleaned values (SGA and z missing whenever GA or BW is missing).
# Parity is top-coded at 8 in the CPP and is capped at 4+ in the analysis,
# so it is unchanged.
#
# Input:  data/analytic_sample_corrected.RDS (07), data/ncppbasa.sas7bdat
# Output: data/analytic_sample_clean.RDS, outputs/implausible_value_counts.csv
# ==============================================================================

library(tidyverse)
library(haven)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
cn <- function(x) { x <- as.character(x); x[x %in% c(".u", ".U", "u", "U")] <- NA; as.numeric(x) }

a <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS"))
stopifnot(nrow(a) == 45268)
raw <- read_sas(file.path(in_dir, "ncppbasa.sas7bdat"),
                col_select = c(MOMID, PREGID, CHILDID, C31, C81, C83, C52)) %>%
  mutate(ht_cm = cn(C81) * 2.54, wt_kg = cn(C83) * 0.453592, cigs = cn(C52), age_raw = cn(C31)) %>%
  select(MOMID, PREGID, CHILDID, ht_cm, wt_kg, cigs, age_raw)
a <- a %>% left_join(raw, by = c("MOMID", "PREGID", "CHILDID"))
stopifnot(nrow(a) == 45268, isTRUE(all.equal(as.numeric(a$age), a$age_raw)))

ex <- !is.na(a$ih1) & a$exclusion == 0
flag <- list(
  `Maternal age under 13 years` = !is.na(a$age) & a$age < 13,
  `Height under 122 cm or over 213 cm` = !is.na(a$ht_cm) & (a$ht_cm < 122 | a$ht_cm > 213),
  `Weight under 30 kg` = !is.na(a$wt_kg) & a$wt_kg < 30,
  `BMI under 13 or over 60` = !is.na(a$bmi) & (a$bmi < 13 | a$bmi > 60),
  `Cigarettes per day coded 70 or 80` = !is.na(a$cigs) & a$cigs %in% c(70, 80),
  `Gestational age outside 20-44 weeks` = !is.na(a$gest_age) & (a$gest_age < 20 | a$gest_age > 44),
  `Birthweight under 300 g or over 6,000 g` = !is.na(a$birthweight) & (a$birthweight < 300 | a$birthweight > 6000))
flag$`Birthweight implausible for gestational age (GA set missing)` <-
  !flag$`Gestational age outside 20-44 weeks` & !flag$`Birthweight under 300 g or over 6,000 g` &
  !is.na(a$gest_age) & !is.na(a$birthweight) &
  ((a$gest_age < 28 & a$birthweight > 2500) | (a$gest_age >= 37 & a$birthweight < 1000))

bmi_bad <- flag$`Height under 122 cm or over 213 cm` | flag$`Weight under 30 kg` | flag$`BMI under 13 or over 60`
counts <- tibble(rule = names(flag), all_45268 = map_int(flag, sum), examined_40700 = map_int(flag, ~ sum(.x & ex))) %>%
  bind_rows(tibble(rule = c("Age set missing (total missing after cleaning)", "BMI set missing (total missing after cleaning)",
                            "Smoking set missing (total missing after cleaning)"),
                   all_45268 = c(sum(is.na(a$age) | flag[[1]]), sum(is.na(a$bmi) | bmi_bad),
                                 sum(is.na(a$smoking) | flag[[5]])),
                   examined_40700 = c(sum((is.na(a$age) | flag[[1]]) & ex), sum((is.na(a$bmi) | bmi_bad) & ex),
                                      sum((is.na(a$smoking) | flag[[5]]) & ex))))
print(counts, width = Inf)

lab <- function(x, l) { attr(x, "label") <- l; x }
clean <- a %>%
  mutate(age     = lab(ifelse(flag[[1]], NA, age), "Maternal age, years (under 13 set missing)"),
         bmi     = lab(ifelse(bmi_bad, NA, bmi), "Pre-pregnancy BMI, kg/m2 (implausible set missing)"),
         smoking = lab(ifelse(flag[[5]], NA, smoking), "0=None, 1=1-19/day, 2=20+/day (codes 70, 80 set missing)"),
         gest_age    = ifelse(flag[[6]] | flag[[8]], NA, gest_age),
         birthweight = ifelse(flag[[7]], NA, birthweight),
         preterm  = as.integer(gest_age < 37),
         preterm3 = case_when(gest_age < 32 ~ 2, gest_age < 37 ~ 1, gest_age >= 37 ~ 0)) %>%
  group_by(gest_age, infant_sex) %>%
  mutate(bw_for_ga_z = (birthweight - mean(birthweight, na.rm = TRUE)) / sd(birthweight, na.rm = TRUE),
         SGA = as.integer(bw_for_ga_z < -1.28)) %>%
  ungroup() %>%
  mutate(bw_for_ga_z = ifelse(is.na(gest_age) | is.na(birthweight), NA, bw_for_ga_z),
         SGA = ifelse(is.na(gest_age) | is.na(birthweight), NA_integer_, SGA)) %>%
  select(-ht_cm, -wt_kg, -cigs, -age_raw)
stopifnot(identical(names(clean), names(readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")))),
          nrow(clean) == 45268)

saveRDS(clean, file.path(in_dir, "analytic_sample_clean.RDS"))
dir.create(out_dir, showWarnings = FALSE)
write_csv(counts, file.path(out_dir, "implausible_value_counts.csv"))
cat("Saved:", file.path(in_dir, "analytic_sample_clean.RDS"), "\n")
