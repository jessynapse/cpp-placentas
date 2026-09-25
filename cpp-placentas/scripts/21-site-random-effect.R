# ==============================================================================
# 21-site-random-effect.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# PART A. Is study site more strongly associated with the exposures or with
#         the outcome? Site variation is put on a common scale with
#         random-intercept logistic models (lme4::glmer), one per variable,
#         in the same 40,700 infants with an observed outcome:
#           y ~ 1 + (1 | site)                 unadjusted
#           y ~ race + age + (1 | site)        adjusted for fully observed
#                                              maternal characteristics
#         for y = any MVM, any AI, definite IH.
#         Reported: site variance, median odds ratio (MOR) and latent-scale
#         ICC, plus the range of prevalence across sites.
#           MOR = exp(sqrt(2 * var) * qnorm(0.75))
#           ICC = var / (var + pi^2 / 3)
#         The MOR is the median odds ratio between two randomly chosen
#         sites for otherwise similar women. MOR = 1 means no site variation.
#
# PART B. Site as a RANDOM effect instead of fixed effects. Mixed Poisson
#         model (log link) with a random intercept for site:
#           ih1 ~ exposure + 11 covariates + (1 | site), weights = IPW
#         fit in each of 20 imputations and pooled with Rubin's rules.
#         Four exposures: MVM any, MVM score, AI any, AI score.
#         Shown next to the fixed-effect primary model (from 14).
#         Notes: (1) model-based SEs from a Poisson model are conservative
#         for a binary outcome. (2) Sibling clustering (MOMID) is not
#         modelled here (few mothers have >1 infant, and GEE results show
#         it barely affects SEs). (3) With 12 sites and large site sizes,
#         random and fixed site effects are expected to give similar
#         exposure estimates. A random intercept also assumes site is
#         unrelated to the exposure, which Part A tests.
#
# Outputs (outputs/):
#   site_variation_exposure_outcome.csv   Part A
#   site_random_effect_models.csv         Part B
# No p-values.
# ==============================================================================

library(tidyverse)
library(haven)
library(lme4)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
n_cores <- min(4, detectCores())   # set to 1 on Windows

# AI score (as 14)
ai_comps <- c("PA02_34", "PA02_35", "PA02_36", "PA02_37",
              "PA02_25", "PA02_26", "PA02_38")
path2 <- read_sas(file.path(in_dir, "path2all.sas7bdat"),
                  col_select = all_of(c("MOMID", "PREGID", "CHILDID", ai_comps))) %>%
  zap_labels() %>%
  mutate(ai_count = rowSums(across(all_of(ai_comps), ~ !is.na(.x) & .x > 0))) %>%
  select(MOMID, PREGID, CHILDID, ai_count)

imp_data <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  filter(.imp > 0, !is.na(ih1), exclusion == 0) %>%
  left_join(path2, by = c("MOMID", "PREGID", "CHILDID")) %>%
  mutate(
    smoking   = factor(smoking,   levels = c("Non-smoker","1-19/day","20+/day")),
    parity    = factor(parity,    levels = c("0","1","2","3","4+")),
    income    = factor(as.character(income),
                       levels = c("$4,000-$5,999","<=$1,999",
                                  "$2,000-$3,999","$6,000-$7,999",
                                  "$8,000-$9,999","$10,000-$14,999",
                                  ">=$15,000")),
    educ      = factor(educ,      levels = c("<HS","Some HS","HS grad","College+")),
    marital   = factor(marital,   levels = c("Married/CL","Single","Wid/Div/Sep")),
    race      = factor(race,      levels = c("White","Black","Puerto Rican","Other")),
    site      = factor(site),
    mvm2_num  = as.numeric(as.character(mvm2)),
    ai_num    = as.numeric(as.character(ai)),
    # centred and scaled continuous covariates help glmer converge; this
    # does not change the exposure coefficients
    age_c     = as.numeric(scale(age)),
    bmi_c     = as.numeric(scale(bmi))
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor)) %>%
  select(.imp, MOMID, ih1, ipw, mvm2_num, mvm, ai_num, ai_count, site, age,
         age_c, bmi_c, race, educ, income, marital, smoking, parity, dm,
         chronic_htn, infant_sex) %>%
  split(.$.imp)
invisible(gc())

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ------------------------------------------------------------------------------
# PART A. SITE VARIATION IN EXPOSURES VS OUTCOME
# Exposures, outcome, race and age are never imputed, so imputation 1 is
# used (identical in all 20).
# ------------------------------------------------------------------------------
d1 <- imp_data[[1]]
vars <- c(mvm2_num = "Any MVM (exposure)", ai_num = "Any AI (exposure)",
          ih1 = "Definite IH (outcome)")

site_var <- function(y, adj) {
  f <- as.formula(paste(y, "~", if (adj) "race + age_c" else "1", "+ (1 | site)"))
  m <- glmer(f, data = d1, family = binomial, control = ctrl)
  v <- as.numeric(VarCorr(m)$site)
  tibble(site_variance = v,
         MOR = exp(sqrt(2 * v) * qnorm(0.75)),
         ICC = v / (v + pi^2 / 3))
}

partA <- imap_dfr(vars, function(lab, y) {
  prev <- d1 %>% group_by(site) %>% summarise(p = 100 * mean(.data[[y]] == 1), .groups = "drop")
  bind_rows(site_var(y, FALSE) %>% mutate(adjustment = "Unadjusted"),
            site_var(y, TRUE)  %>% mutate(adjustment = "Adjusted for race and age")) %>%
    mutate(variable = lab,
           overall_prevalence = sprintf("%.1f%%", 100 * mean(d1[[y]] == 1)),
           site_range = sprintf("%.1f%% to %.1f%%", min(prev$p), max(prev$p)),
           max_min_ratio = round(max(prev$p) / min(prev$p), 1))
}) %>%
  mutate(across(c(site_variance, MOR), ~ round(.x, 2)), ICC = round(ICC, 3)) %>%
  select(variable, adjustment, overall_prevalence, site_range, max_min_ratio,
         site_variance, MOR, ICC)

cat("\n==================== PART A: SITE VARIATION ====================\n")
print(partA, width = Inf)

# ------------------------------------------------------------------------------
# PART B. SITE AS A RANDOM EFFECT
# ------------------------------------------------------------------------------
confounders <- "age_c + bmi_c + race + educ + income + marital +
                smoking + parity + dm + chronic_htn + infant_sex"
exposures <- c(mvm2_num = "MVM any vs none", mvm = "MVM score, per 1 point",
               ai_num = "AI any vs none", ai_count = "AI score, per 1 compartment")

fit_re <- function(e) {
  f <- as.formula(paste("ih1 ~", e, "+", confounders, "+ (1 | site)"))
  one <- function(d) tryCatch({
    m <- suppressWarnings(glmer(f, data = d, family = poisson(link = "log"),
                                weights = ipw, control = ctrl))
    list(b = fixef(m)[[e]], v = vcov(m)[e, e],
         s2 = as.numeric(VarCorr(m)$site))
  }, error = function(err) NULL)
  fits <- mclapply(imp_data, one, mc.cores = n_cores)
  redo <- which(sapply(fits, is.null))
  for (i in redo) fits[[i]] <- one(imp_data[[i]])
  fits <- Filter(Negate(is.null), fits)
  b <- sapply(fits, `[[`, "b"); v <- sapply(fits, `[[`, "v")
  m <- length(b)
  se <- sqrt(mean(v) + (1 + 1/m) * var(b))
  tibble(exposure = exposures[[e]], rr = exp(mean(b)),
         lci = exp(mean(b) - 1.96 * se), uci = exp(mean(b) + 1.96 * se),
         site_variance_ih = mean(sapply(fits, `[[`, "s2")), n_imps = m)
}

partB <- map_dfr(names(exposures), function(e) { cat(" ", exposures[[e]], "\n"); fit_re(e) }) %>%
  mutate(`Site random effect + IPW` = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci),
         MRR_site = round(exp(sqrt(2 * site_variance_ih) * qnorm(0.75)), 2))

fixed <- read_csv(file.path(out_dir, "four_exposures_results.csv"), show_col_types = FALSE) %>%
  mutate(exposure = recode(exposure,
                           "MVM score, per 1 point (0 to 8)" = "MVM score, per 1 point",
                           "AI score, per 1 compartment (0 to 7)" = "AI score, per 1 compartment")) %>%
  select(exposure, model, rr_ci) %>%
  pivot_wider(names_from = model, values_from = rr_ci)

partB_out <- partB %>%
  left_join(fixed, by = "exposure") %>%
  transmute(exposure,
            n_cases = sum(d1$ih1 == 1), n_obs = nrow(d1),
            Crude,
            `Primary: site as fixed effect + IPW` = `Primary: + site + IPW`,
            `Site as random effect + IPW` = `Site random effect + IPW`,
            `Secondary: no site + IPW` = `Secondary: no site + IPW`,
            median_rate_ratio_site = MRR_site)

cat("\n==================== PART B: SITE AS RANDOM EFFECT ====================\n")
print(partB_out, width = Inf)

write_csv(partA,     file.path(out_dir, "site_variation_exposure_outcome.csv"))
write_csv(partB_out, file.path(out_dir, "site_random_effect_models.csv"))
cat("\nSaved: outputs/site_variation_exposure_outcome.csv, site_random_effect_models.csv\n")
