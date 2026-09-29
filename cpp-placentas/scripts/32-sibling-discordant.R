# ==============================================================================
# 32-sibling-discordant.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Feasibility of a sibling (within-mother) analysis, suggested by Alexa
# Freedman ("Probably not feasible given low prevalence of IH, but could
# look at discordant siblings to improve control for maternal
# characteristics?"). Supplemental.
#
# A within-mother comparison controls for everything shared by a mother's
# pregnancies (genetics, stable socioeconomic and lifestyle factors, race,
# and study site, which never changes within a mother in CPP). Only
# mothers whose pregnancies differ on BOTH IH and the exposure contribute
# information to the estimate.
#
# Part 1. Feasibility counts in the 40,700 infants with an observed outcome
# Part 2. Conditional logistic regression (survival::clogit), stratified on
#         MOMID, for any MVM and any AI:
#           crude
#           + pregnancy-varying covariates: maternal age, parity, smoking,
#             infant sex (complete case on these, observed data only)
#         For comparison, a standard (between-mother) logistic GEE in the
#         SAME sibling sample, to separate the effect of restricting the
#         sample from the effect of the within-mother comparison.
#         With IH at 1.7%, the odds ratio approximates the risk ratio.
#
# Outputs (outputs/):
#   sibling_feasibility.csv   counts
#   sibling_models.csv        within- vs between-mother estimates
# No p-values.
# ==============================================================================

library(tidyverse)
library(survival)
library(geepack)

in_dir  <- "/Users/wongjj/Downloads"   # folder with analytic_sample_corrected.RDS
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE)

fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  filter(!is.na(ih1), exclusion == 0) %>%
  mutate(mvm_any = as.integer(as.character(mvm2) == "1"),
         ai_any  = as.integer(as.character(ai) == "1"),
         smoking = factor(smoking), parity = factor(parity),
         infant_sex = factor(infant_sex))
stopifnot(nrow(d) == 40700, sum(d$ih1) == 700)
stopifnot(all(tapply(d$site, d$MOMID, n_distinct) == 1))

# ------------------------------------------------------------------------------
# PART 1. FEASIBILITY
# ------------------------------------------------------------------------------
sib <- d %>% group_by(MOMID) %>% filter(n() > 1) %>% ungroup()
disc_ih <- sib %>% group_by(MOMID) %>% filter(n_distinct(ih1) > 1) %>% ungroup()
informative <- function(x) disc_ih %>% group_by(MOMID) %>%
  filter(n_distinct(.data[[x]]) > 1) %>% ungroup()
inf_mvm <- informative("mvm_any"); inf_ai <- informative("ai_any")

row <- function(label, dd) tibble(
  step = label, mothers = fmt_n(n_distinct(dd$MOMID)),
  pregnancies = fmt_n(nrow(dd)), `IH cases` = fmt_n(sum(dd$ih1)))

feas <- bind_rows(
  row("All infants with observed IH", d),
  row("Mothers with 2 or more infants in the sample", sib),
  row("  IH-discordant mothers (at least one case and one non-case)", disc_ih),
  row("    also discordant on any MVM (informative for MVM)", inf_mvm),
  row("    also discordant on any AI (informative for AI)", inf_ai))

# Within informative mothers: is IH in the exposed or unexposed pregnancy?
dir_tab <- function(dd, x) dd %>% filter(ih1 == 1) %>%
  summarise(exposed = sum(.data[[x]] == 1), unexposed = sum(.data[[x]] == 0)) %>%
  mutate(exposure = x)
direction <- bind_rows(dir_tab(inf_mvm, "mvm_any"), dir_tab(inf_ai, "ai_any"))

# ------------------------------------------------------------------------------
# PART 2. WITHIN-MOTHER (clogit) VS BETWEEN-MOTHER (GEE) IN SIBLING SAMPLE
# ------------------------------------------------------------------------------
covs <- "age + parity + smoking + infant_sex"
sib_cc <- sib %>% filter(if_all(c(age, parity, smoking, infant_sex), ~ !is.na(.x))) %>%
  arrange(MOMID)

est <- function(fit, term) {
  b <- coef(fit)[term]; se <- sqrt(diag(vcov(fit)))[term]
  sprintf("%.2f (%.2f, %.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
}

models <- imap_dfr(c(mvm_any = "Any MVM", ai_any = "Any AI"), function(lab, x) {
    fits <- list(
      `Within-mother (conditional logistic), crude` =
        clogit(as.formula(paste("ih1 ~", x, "+ strata(MOMID)")), data = sib_cc),
      `Within-mother (conditional logistic), + age, parity, smoking, sex` =
        clogit(as.formula(paste("ih1 ~", x, "+", covs, "+ strata(MOMID)")), data = sib_cc),
      `Between-mother (logistic GEE), crude, same sibling sample` =
        geeglm(as.formula(paste("ih1 ~", x)), data = sib_cc, family = binomial,
               id = MOMID, corstr = "exchangeable"),
      `Between-mother (logistic GEE), + age, parity, smoking, sex, same sibling sample` =
        geeglm(as.formula(paste("ih1 ~", x, "+", covs)), data = sib_cc,
               family = binomial, id = MOMID, corstr = "exchangeable"))
    imap_dfr(fits, ~ tibble(exposure = lab, model = .y, `OR (95% CI)` = est(.x, x)))
  }) %>%
  mutate(n_cases = sum(sib_cc$ih1), n_obs = nrow(sib_cc),
         n_mothers = n_distinct(sib_cc$MOMID))

cat("\n==================== SIBLING FEASIBILITY ====================\n")
print(feas, width = Inf)
cat("\nIH cases in informative mothers, by exposure of the affected sibling:\n")
print(direction)
cat("\n==================== WITHIN- VS BETWEEN-MOTHER (definite IH) ====================\n")
print(models, width = Inf)

write_csv(feas, file.path(out_dir, "sibling_feasibility.csv"))
write_csv(models, file.path(out_dir, "sibling_models.csv"))
cat("\nSaved: outputs/sibling_feasibility.csv, outputs/sibling_models.csv\n")
