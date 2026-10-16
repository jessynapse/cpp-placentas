# ==============================================================================
# 38-included-excluded.R
# CPP Placental Pathology -> Infantile Hemangioma
# October 2026
#
# Selection bias table required by Paediatric and Perinatal Epidemiology:
# characteristics of the 40,700 infants examined at one year (included)
# compared with the 4,568 excluded because they died before the one-year
# examination (1,834) or survived without an IH examination (2,734).
# All 45,268 pregnancies in the analytic sample, observed data.
# Same variable list and format as Table 1 (33), plus placental pathology
# and study site.
#
# Output: outputs/etable_included_excluded.csv
# ==============================================================================

library(tidyverse)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
n_ <- function(x) format(x, big.mark = ",", trim = TRUE)
pct_n <- function(k, n) sprintf("%.1f%% (%s)", 100 * k / n, n_(k))

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS"))
site_names <- c("5" = "Boston", "10" = "Buffalo", "15" = "New Orleans",
                "31" = "New York (Columbia)", "37" = "Baltimore", "45" = "Virginia",
                "50" = "Minnesota", "55" = "New York Medical College", "60" = "Oregon",
                "66" = "Pennsylvania", "71" = "Providence", "82" = "Tennessee")

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


d <- d %>% mutate(site_lab = factor(paste(as.character(site), site_names[as.character(site)]),
                                    levels = paste(names(site_names), site_names)))
vars <- c(list(list("Placental pathology", NA, "header"),
               list("Maternal vascular malperfusion", "mvm_any", "cat"),
               list("Acute inflammation", "ai_any", "cat"),
               list("Joint MVM and AI", "joint", "cat"),
               list("Study site", "site_lab", "cat"),
               list("Maternal and infant characteristics", NA, "header")), vars)

grp <- list(
  `Included (examined at one year)` = filter(d, !is.na(ih1), exclusion == 0),
  `Excluded, all` = filter(d, is.na(ih1) | exclusion == 1),
  `Excluded, died before one year` = filter(d, exclusion == 1),
  `Excluded, no one-year examination` = filter(d, exclusion == 0, is.na(ih1)))
stopifnot(nrow(grp[[1]]) == 40700, nrow(grp[[2]]) == 4568)
tab <- build_table(grp)
tab <- bind_rows(tibble(Characteristic = "", !!!as.list(map_chr(grp, ~ sprintf("n = %s", n_(nrow(.x)))))), tab)
print(tab, n = Inf, width = Inf)
write_csv(tab, file.path(out_dir, "etable_included_excluded.csv"))
cat("\nSaved: outputs/etable_included_excluded.csv\n")
