# ==============================================================================
# 22-site-coefficients.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# The study site terms from the PRIMARY adjusted model: adjusted RR of
# definite IH for each site compared with Boston (site 5, the largest
# site), holding the exposure, 11 covariates and IPW constant.
# Shown for the MVM any and AI any primary models.
#
# Same methods as the primary model (10/14): 40,700 infants, 700 cases,
# 20 imputations, Rubin's rules, modified Poisson GEE (exchangeable,
# clustered on MOMID), IPW from 09.
#
# Output: outputs/site_coefficients.csv. No p-values.
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
n_cores <- min(4, detectCores())   # set to 1 on Windows

site_levels <- c("5", "10", "15", "31", "37", "45", "50", "55", "60", "66", "71", "82")

imp_data <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  filter(.imp > 0, !is.na(ih1), exclusion == 0) %>%
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
    site      = factor(as.character(site), levels = site_levels),  # Boston = reference
    mvm2      = factor(mvm2),
    ai        = factor(ai)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor)) %>%
  select(.imp, MOMID, ih1, ipw, mvm2, ai, site, age, bmi, race, educ, income,
         marital, smoking, parity, dm, chronic_htn, infant_sex) %>%
  arrange(MOMID) %>%
  split(.$.imp)
invisible(gc())

pool_gee <- function(fits) {
  m        <- length(fits)
  coef_mat <- do.call(rbind, lapply(fits, coef))
  se_mat   <- do.call(rbind, lapply(fits, function(f) sqrt(diag(vcov(f)))))
  Q_bar <- colMeans(coef_mat)
  SE    <- sqrt(colMeans(se_mat^2) + (1 + 1/m) * apply(coef_mat, 2, var))
  tibble(term = names(Q_bar), rr = exp(Q_bar),
         lci = exp(Q_bar - 1.96 * SE), uci = exp(Q_bar + 1.96 * SE))
}
confounders <- "age + bmi + race + educ + income + marital +
                smoking + parity + dm + chronic_htn + infant_sex"

run <- function(e) {
  f <- as.formula(paste("ih1 ~", e, "+", confounders, "+ site"))
  fit_one <- function(d) tryCatch(
    geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
           corstr = "exchangeable", weights = ipw),
    error = function(err) err)
  fits <- mclapply(imp_data, fit_one, mc.cores = n_cores)
  redo <- which(!sapply(fits, function(x) inherits(x, "geeglm")))
  for (i in redo) fits[[i]] <- fit_one(imp_data[[i]])
  fits <- fits[sapply(fits, function(x) inherits(x, "geeglm"))]
  pool_gee(fits) %>% filter(grepl("^site", term)) %>%
    mutate(site = sub("^site", "", term),
           rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci))
}

d1 <- imp_data[[1]]
counts <- d1 %>% group_by(site) %>%
  summarise(n_cases = sum(ih1 == 1), n_obs = n(), .groups = "drop") %>%
  mutate(site = as.character(site),
         ih_pct = sprintf("%.1f%% (%s/%s)", 100 * n_cases / n_obs,
                          format(n_cases, big.mark = ",", trim = TRUE),
                          format(n_obs, big.mark = ",", trim = TRUE)))

res <- bind_rows(run("mvm2") %>% mutate(model = "Primary model, MVM any"),
                 run("ai")   %>% mutate(model = "Primary model, AI any"))

out <- counts %>%
  left_join(res %>% select(site, model, rr_ci) %>%
              pivot_wider(names_from = model, values_from = rr_ci), by = "site") %>%
  mutate(across(starts_with("Primary"), ~ ifelse(site == "5", "1.00 (ref)", .x)),
         site = factor(site, site_levels)) %>%
  arrange(site) %>%
  mutate(site = ifelse(site == "5", "5 (Boston, reference)", as.character(site)))

cat("\n==================== SITE TERMS IN THE PRIMARY MODEL (definite IH) ====================\n")
print(out, n = Inf, width = Inf)
write_csv(out, file.path(out_dir, "site_coefficients.csv"))
cat("\nSaved: outputs/site_coefficients.csv\n")
