# ==============================================================================
# 26-site-global-association.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Overall association of study site AS ONE VARIABLE (all 12 levels) with
# any MVM, any AI and definite IH, in the 40,700 infants with an observed
# outcome. Requested by Stefanie Hinkle: is site more strongly associated
# with the placental exposures or with IH?
#
# For each of MVM ~ site, AI ~ site, IH ~ site (logistic regression):
#   - Likelihood ratio chi-square for the site variable (11 df):
#       unadjusted:              y ~ site         vs  y ~ 1
#       adjusted (race, age):    y ~ race + age + site  vs  y ~ race + age
#   - Cramér's V for the site x variable table (0 = none, 1 = perfect)
#   - McFadden pseudo R-squared for site alone
#
# P-values are not reported (CLAUDE.md). With 11 df, a chi-square above
# 19.7 corresponds to the conventional 0.05 threshold.
#
# Output: outputs/site_global_association.csv
# ==============================================================================

library(tidyverse)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  filter(!is.na(ih1), exclusion == 0) %>%
  mutate(site = factor(site), race = factor(race),
         mvm_any = as.integer(mvm2 == 1), ai_any = as.integer(ai == 1),
         ih = as.integer(ih1 == 1))
stopifnot(nrow(d) == 40700)

vars <- c(mvm_any = "Any MVM (placenta)", ai_any = "Any AI (placenta)",
          ih = "Definite IH (outcome)")

res <- imap_dfr(vars, function(lab, y) {
  fit <- function(rhs) glm(as.formula(paste(y, "~", rhs)), data = d, family = binomial)
  m0 <- fit("1"); ms <- fit("site"); ma <- fit("race + age"); mas <- fit("race + age + site")
  lr_unadj <- deviance(m0) - deviance(ms)
  lr_adj   <- deviance(ma) - deviance(mas)
  tab <- table(d$site, d[[y]])
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  V <- sqrt(chi / (sum(tab) * (min(dim(tab)) - 1)))
  tibble(variable = lab,
         prevalence = sprintf("%.1f%%", 100 * mean(d[[y]])),
         `LR chi-square for site, unadjusted (11 df)` = round(lr_unadj, 1),
         `LR chi-square for site, adjusted for race and age (11 df)` = round(lr_adj, 1),
         `Cramer's V` = round(as.numeric(V), 3),
         `McFadden pseudo R2, site alone` = round(1 - deviance(ms) / deviance(m0), 3))
})

cat("\n==================== SITE AS ONE VARIABLE ====================\n")
print(res, width = Inf)
write_csv(res, file.path(out_dir, "site_global_association.csv"))
cat("\nSaved: outputs/site_global_association.csv\n")
