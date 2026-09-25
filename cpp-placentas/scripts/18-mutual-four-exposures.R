# ==============================================================================
# 18-mutual-four-exposures.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Mutually adjusted MVM and AI models for all four exposures:
#   Binary pair:      MVM any + AI any in the same model
#   Continuous pair:  MVM score (per 1 point) + AI score (per 1 compartment)
#                     in the same model
# Each exposure is adjusted for its counterpart of the same type. All four
# are not put in one model, because MVM any and MVM score (and AI any and
# AI score) measure the same exposure.
#
# The four-level joint exposure (neither, MVM only, AI only, both) is binary
# only and is in 16-joint-mvm-ai.R.
#
# Same methods as 13-17: definite IH, 40,700 infants, 700 cases, 20
# imputations, Rubin's rules, modified Poisson GEE (exchangeable, MOMID).
# Models: crude (both exposures, no covariates), secondary (11 covariates
# + IPW), primary (+ site + IPW). The single-exposure primary estimates
# from 14 are shown alongside for comparison.
#
# Output: outputs/mutual_4exp.csv. No p-values.
# ==============================================================================

library(tidyverse)
library(haven)
library(geepack)
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
    mvm2      = factor(mvm2),
    ai        = factor(ai)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor)) %>%
  select(.imp, MOMID, ih1, ipw, mvm2, mvm, ai, ai_count, site, age, bmi, race,
         educ, income, marital, smoking, parity, dm, chronic_htn, infant_sex) %>%
  arrange(MOMID) %>%
  split(.$.imp)
stopifnot(!any(is.na(imp_data[[1]]$ai_count)))
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
models <- c("Crude, mutually adjusted" = NA,
            "Secondary: + covariates + IPW, no site" = "",
            "Primary: + covariates + site + IPW" = "+ site")
pairs  <- c("Binary pair" = "mvm2 + ai", "Continuous pair" = "mvm + ai_count")
labels <- c(mvm21 = "MVM any vs none", mvm = "MVM score, per 1 point",
            ai1 = "AI any vs none", ai_count = "AI score, per 1 compartment")

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
  pool_gee(fits[!bad]) %>% filter(term %in% names(labels))
}

res <- imap_dfr(pairs, function(expo, pr) imap_dfr(models, function(extra, lab) {
  cat(" ", pr, "/", lab, "\n")
  run(expo, extra) %>% mutate(pair = pr, model = lab)
})) %>%
  mutate(exposure = factor(labels[term], levels = labels),
         rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci))

# single-exposure primary estimates from 14, for comparison
single <- read_csv(file.path(out_dir, "four_exposures_results.csv"), show_col_types = FALSE) %>%
  filter(model == "Primary: + site + IPW") %>%
  mutate(exposure = factor(labels[c("MVM any vs none" = "mvm21",
                                    "MVM score, per 1 point (0 to 8)" = "mvm",
                                    "AI any vs none" = "ai1",
                                    "AI score, per 1 compartment (0 to 7)" = "ai_count")[exposure]],
                           levels = labels)) %>%
  transmute(exposure, `Primary, single exposure (for comparison)` = rr_ci)

out <- res %>%
  select(exposure, model, rr_ci) %>%
  pivot_wider(names_from = model, values_from = rr_ci) %>%
  left_join(single, by = "exposure") %>%
  arrange(exposure) %>%
  mutate(adjusted_for = recode(as.character(exposure),
                               "MVM any vs none" = "AI any",
                               "AI any vs none" = "MVM any",
                               "MVM score, per 1 point" = "AI score",
                               "AI score, per 1 compartment" = "MVM score"),
         n_cases = sum(imp_data[[1]]$ih1 == 1), n_obs = nrow(imp_data[[1]])) %>%
  relocate(adjusted_for, .after = exposure)

cat("\n==================== MUTUALLY ADJUSTED, FOUR EXPOSURES ====================\n")
print(out, width = Inf)
write_csv(out, file.path(out_dir, "mutual_4exp.csv"))
cat("\nSaved: outputs/mutual_4exp.csv\n")
