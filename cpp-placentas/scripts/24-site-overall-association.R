# ==============================================================================
# 24-site-overall-association.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Overall strength of the association between study site (all 12 levels)
# and each of: any MVM, any AI (placenta) and definite IH (outcome), on the
# relative risk scale, to answer "is site more strongly associated with
# the placental exposures or with IH?" (Stefanie Hinkle).
#
# For each variable, in the 40,700 infants with an observed outcome:
#   1. Range of site RRs vs Boston, adjusted for race and age (from 23)
#   2. Median rate ratio (MRR) from a random-intercept log-link (Poisson)
#      model, y ~ (1 | site), unadjusted and + race + age:
#        MRR = exp(sqrt(2 * var) * qnorm(0.75))
#      the median RR between two randomly chosen sites (1 = no variation)
#   3. % of deviance explained by site in log-link models:
#        site alone:              (D_null - D_site) / D_null
#        beyond race and age:     (D_race_age - D_race_age_site) / D_null
#      Poisson log-link deviance is used for all three variables so the
#      measure is on the same scale (log-binomial does not converge for
#      MVM).
#
# Output: outputs/site_overall_association.csv. No p-values.
# ==============================================================================

library(tidyverse)
library(lme4)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")

site_levels <- c("5", "10", "15", "31", "37", "45", "50", "55", "60", "66", "71", "82")
d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  filter(!is.na(ih1), exclusion == 0) %>%
  mutate(site = factor(as.character(site), levels = site_levels),
         race = factor(race, levels = c(1, 2, 4, 8)),
         age_c = (age - mean(age)) / sd(age),
         mvm_any = as.integer(mvm2 == 1),
         ai_any  = as.integer(ai == 1),
         ih      = as.integer(ih1 == 1))
stopifnot(nrow(d) == 40700)

vars <- c(mvm_any = "Any MVM (placenta)", ai_any = "Any AI (placenta)",
          ih = "Definite IH (outcome)")
ctrl <- glmerControl(optimizer = "bobyqa")

# site RR ranges, adjusted for race and age (from 23)
rng <- read_csv(file.path(out_dir, "site_log_binomial.csv"), show_col_types = FALSE)
rr_range <- function(col) {
  x <- as.numeric(sub(" .*", "", rng[[col]][rng$site != "5 (Boston, ref)"]))
  sprintf("%.2f to %.2f", min(c(x, 1)), max(c(x, 1)))
}
range_cols <- c(mvm_any = "Any MVM RR, adjusted for race and age",
                ai_any  = "Any AI RR, adjusted for race and age",
                ih      = "Definite IH RR, adjusted for race and age")

res <- imap_dfr(vars, function(lab, y) {
  mrr <- function(adj) {
    f <- as.formula(paste(y, "~", if (adj) "race + age_c" else "1", "+ (1 | site)"))
    v <- as.numeric(VarCorr(glmer(f, data = d, family = poisson(link = "log"),
                                  control = ctrl))$site)
    exp(sqrt(2 * v) * qnorm(0.75))
  }
  dev <- function(rhs) deviance(glm(as.formula(paste(y, "~", rhs)), data = d,
                                    family = poisson(link = "log")))
  D0 <- dev("1"); Ds <- dev("site"); Dra <- dev("race + age_c"); Dras <- dev("race + age_c + site")
  tibble(variable = lab,
         prevalence = sprintf("%.1f%%", 100 * mean(d[[y]])),
         `Site RR range vs Boston (adj. race, age)` = rr_range(range_cols[[y]]),
         `MRR, unadjusted` = round(mrr(FALSE), 2),
         `MRR, adjusted for race and age` = round(mrr(TRUE), 2),
         `% deviance explained by site alone` = round(100 * (D0 - Ds) / D0, 1),
         `% deviance explained by site beyond race and age` = round(100 * (Dra - Dras) / D0, 1))
})

cat("\n==================== OVERALL ASSOCIATION OF SITE (12 levels) ====================\n")
print(res, width = Inf)
write_csv(res, file.path(out_dir, "site_overall_association.csv"))
cat("\nSaved: outputs/site_overall_association.csv\n")
