# ==============================================================================
# 25-site-logistic.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Site-specific ODDS RATIOS for the exposures and for the outcome, from
# logistic regression (Boston = reference), as requested by Stefanie Hinkle.
#   any MVM ~ site        any AI ~ site        definite IH ~ site
# each unadjusted and adjusted for race and maternal age (both fully
# observed, so no imputation is needed).
#
# Sample: the 40,700 infants with an observed IH outcome (same people in
# every model). Exposures and outcome are never imputed.
#
# Logistic regression with robust (sandwich) SEs clustered on MOMID
# (geepack, independence working correlation) to account for siblings.
# Companion to 23, which gives relative risks from log-binomial models.
#
# Output: outputs/site_logistic.csv. No p-values.
# ==============================================================================

library(tidyverse)
library(geepack)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")

site_levels <- c("5", "10", "15", "31", "37", "45", "50", "55", "60", "66", "71", "82")

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  filter(!is.na(ih1), exclusion == 0) %>%
  mutate(site = factor(as.character(site), levels = site_levels),
         race = factor(race, levels = c(1, 2, 4, 8),
                       labels = c("White", "Black", "Puerto Rican", "Other")),
         age_c = age - mean(age),
         mvm_any = as.integer(mvm2 == 1),
         ai_any  = as.integer(ai == 1)) %>%
  arrange(MOMID)
stopifnot(nrow(d) == 40700, !any(is.na(d$site)), !any(is.na(d$race)), !any(is.na(d$age)))

fmt <- function(x, n) sprintf("%.1f%% (%s/%s)", 100 * x / n,
                              format(x, big.mark = ",", trim = TRUE),
                              format(n, big.mark = ",", trim = TRUE))

fit_rr <- function(y, adj) {
  f <- as.formula(paste(y, "~ site", if (adj) "+ race + age_c" else ""))
  m <- geeglm(f, data = d, family = binomial(link = "logit"), id = MOMID,
              corstr = "independence")
  b <- coef(m); se <- sqrt(diag(vcov(m)))
  tibble(term = names(b), or = exp(b), lci = exp(b - 1.96 * se), uci = exp(b + 1.96 * se)) %>%
    filter(grepl("^site", term)) %>%
    transmute(site = sub("^site", "", term),
              rr_ci = sprintf("%.2f (%.2f, %.2f)", or, lci, uci),
              method = "logistic")
}

vars <- c(mvm_any = "Any MVM", ai_any = "Any AI", ih1 = "Definite IH")

prev <- d %>% group_by(site) %>%
  summarise(n = n(),
            `Any MVM %` = fmt(sum(mvm_any), n()),
            `Any AI %`  = fmt(sum(ai_any), n()),
            `Definite IH %` = fmt(sum(ih1 == 1), n()),
            .groups = "drop") %>%
  mutate(site = as.character(site))

rr_tab <- imap_dfr(vars, function(lab, y) bind_rows(
  fit_rr(y, FALSE) %>% mutate(var = lab, adj = "unadjusted"),
  fit_rr(y, TRUE)  %>% mutate(var = lab, adj = "adjusted for race and age")))


out <- prev %>%
  left_join(rr_tab %>% mutate(col = paste0(var, " OR, ", adj)) %>%
              select(site, col, rr_ci) %>%
              pivot_wider(names_from = col, values_from = rr_ci), by = "site") %>%
  mutate(across(contains(" OR, "), ~ ifelse(site == "5", "1.00 (ref)", .x)),
         site = factor(site, site_levels)) %>%
  arrange(site) %>%
  mutate(site = ifelse(site == "5", "5 (Boston, ref)", as.character(site))) %>%
  select(site, n,
         `Any MVM %`, `Any MVM OR, unadjusted`, `Any MVM OR, adjusted for race and age`,
         `Any AI %`, `Any AI OR, unadjusted`, `Any AI OR, adjusted for race and age`,
         `Definite IH %`, `Definite IH OR, unadjusted`, `Definite IH OR, adjusted for race and age`)

cat("\n==================== SITE ODDS RATIOS (logistic) (Boston = reference) ====================\n")
print(out, n = Inf, width = Inf)
write_csv(out, file.path(out_dir, "site_logistic.csv"))
cat("\nSaved: outputs/site_logistic.csv\n")
