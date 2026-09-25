# ==============================================================================
# 20-grade-stage-separate.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Dose-response for MVM and AI SEPARATELY, each in 3 categories with
# none as the reference (each exposure in its own model):
#   MVM grade: none (0), low (score 1-2), high (score 3+)   [mvm3]
#   AI stage:  none, low, high (high = amnion or umbilical
#              artery involved, per Freedman et al. 2024)  [ai3]
#
#   1. MVM grade: none (ref), low, high
#   2. AI stage:  none (ref), low, high
# Weights are the primary weights from 09 (denominator includes mvm2 and
# ai), so estimates are comparable with the other tables. 12 rebuilt the
# weights with ai3 instead; results differ by at most 0.01.
#
# Same methods as 13-18: definite IH, 40,700 infants, 700 cases, 20
# imputations, Rubin's rules, modified Poisson GEE (exchangeable, MOMID).
# Models: crude, secondary (11 covariates + IPW), primary (+ site + IPW).
#
# Output: outputs/grade_stage_separate.csv
# No p-values.
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
n_cores <- min(4, detectCores())   # set to 1 on Windows

fmt <- function(x, n) sprintf("%.1f%% (%s/%s)", 100 * x / n,
                              format(x, big.mark = ",", trim = TRUE),
                              format(n, big.mark = ",", trim = TRUE))
lv <- c("None", "Low", "High")

imp_data <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  filter(.imp > 0, !is.na(ih1), exclusion == 0) %>%
  mutate(
    mvm3 = factor(as.character(mvm3), levels = lv),
    ai3  = factor(as.character(ai3),  levels = lv),
    joint = factor(paste0("MVM ", mvm3, " / AI ", ai3),
                   levels = as.vector(outer(lv, lv, function(m, a) paste0("MVM ", m, " / AI ", a)))),
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
    site      = factor(site)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor)) %>%
  select(.imp, MOMID, ih1, ipw, mvm3, ai3, joint, site, age, bmi, race, educ,
         income, marital, smoking, parity, dm, chronic_htn, infant_sex) %>%
  arrange(MOMID) %>%
  split(.$.imp)
stopifnot(!any(is.na(imp_data[[1]]$joint)))
invisible(gc())

# ------------------------------------------------------------------------------
# COUNTS
# ------------------------------------------------------------------------------
d1 <- imp_data[[1]]

# ------------------------------------------------------------------------------
# MODELS
# ------------------------------------------------------------------------------
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
models <- c("Crude" = NA, "Secondary: covariates + IPW, no site" = "",
            "Primary: covariates + site + IPW" = "+ site")

run <- function(expo, extra) {
  rhs <- if (is.na(extra)) expo else paste(expo, "+", confounders, extra)
  f <- as.formula(paste("ih1 ~", rhs))
  use_ipw <- !is.na(extra)
  fit_one <- function(d) {
    d$.w <- if (use_ipw) d$ipw else rep(1, nrow(d))
    tryCatch(geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                    corstr = "exchangeable", weights = .w),
             error = function(err) err)
  }
  fits <- mclapply(imp_data, fit_one, mc.cores = n_cores)
  redo <- which(!sapply(fits, function(x) inherits(x, "geeglm")))
  for (i in redo) fits[[i]] <- fit_one(imp_data[[i]])
  bad <- sapply(fits, function(x) !inherits(x, "geeglm"))
  if (any(bad)) message("  ", sum(bad), " fit(s) failed: ", rhs)
  pool_gee(fits[!bad]) %>% filter(grepl("^joint|^mvm3|^ai3", term))
}

specs <- c("MVM grade" = "mvm3", "AI stage" = "ai3")
res <- imap_dfr(specs, function(expo, sp) imap_dfr(models, function(extra, lab) {
  cat(" ", sp, "/", lab, "\n")
  run(expo, extra) %>% mutate(analysis = sp, model = lab)
})) %>%
  mutate(level = sub("^mvm3|^ai3", "", term),
         rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci))

lvl_counts <- bind_rows(
  d1 %>% group_by(level = as.character(mvm3)) %>%
    summarise(n_cases = sum(ih1 == 1), n_obs = n(), .groups = "drop") %>%
    mutate(analysis = "MVM grade"),
  d1 %>% group_by(level = as.character(ai3)) %>%
    summarise(n_cases = sum(ih1 == 1), n_obs = n(), .groups = "drop") %>%
    mutate(analysis = "AI stage")) %>%
  mutate(ih_pct = fmt(n_cases, n_obs))

out <- lvl_counts %>%
  left_join(res %>% select(analysis, level, model, rr_ci) %>%
              pivot_wider(names_from = model, values_from = rr_ci),
            by = c("analysis", "level")) %>%
  mutate(across(all_of(names(models)), ~ ifelse(level == "None", "1.00 (ref)", .x)),
         analysis = factor(analysis, names(specs)), level = factor(level, lv)) %>%
  arrange(analysis, level) %>%
  select(analysis, level, ih_pct, n_cases, n_obs, `Crude`,
         `Primary: covariates + site + IPW`, `Secondary: covariates + IPW, no site`)

cat("\n==================== MVM GRADE AND AI STAGE, SEPARATE MODELS ====================\n")
print(out %>% select(-n_cases, -n_obs), n = Inf, width = Inf)
write_csv(out, file.path(out_dir, "grade_stage_separate.csv"))
cat("\nSaved: outputs/grade_stage_separate.csv\n")
