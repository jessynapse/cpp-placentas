# ==============================================================================
# 46-verify-any-mvm-any-ai.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Independent check of the primary-model estimates for any MVM and any AI,
# which both round to aRR 1.14 (0.96, 1.35) in 37 (Ellen Francis queried
# whether identical estimates indicate an error).
#
# Rebuilt from the raw data (07 -> 35 -> 36), then re-fitted here with
# separately written code (not the fit_pool() function in 37):
#   1. Confirms MVM and AI are different variables (cross-tab, agreement)
#   2. Primary model for each exposure in each of the 50 imputations,
#      pooled with Rubin's rules, reported to four decimals
#   3. Per-imputation estimates, to show they differ across imputations
#   4. The same models without study site (secondary), without weights,
#      and crude, to show the two exposures diverge once anything changes
#
# Input:  data/weighted_long_m50.RDS (36)
# Output: outputs/verify_any_mvm_ai.csv, outputs/verify_any_mvm_ai_by_imp.csv
# ==============================================================================

library(tidyverse)
library(geepack)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
if (basename(getwd()) == "data") out_dir <- "../outputs"

w <- readRDS(file.path(in_dir, "weighted_long_m50.RDS"))

d <- w %>%
  filter(.imp > 0, !is.na(ih1), exclusion == 0) %>%
  mutate(
    mvm_any = as.numeric(as.character(mvm2)),
    ai_any  = as.numeric(as.character(ai)),
    smoking = factor(as.character(smoking), levels = c("Non-smoker", "1-19/day", "20+/day")),
    parity  = factor(as.character(parity),  levels = c("0", "1", "2", "3", "4+")),
    income  = factor(as.character(income),
                     levels = c("$4,000-$5,999", "<=$1,999", "$2,000-$3,999", "$6,000-$7,999",
                                "$8,000-$9,999", "$10,000-$14,999", ">=$15,000")),
    educ    = factor(as.character(educ),    levels = c("<HS", "Some HS", "HS grad", "College+")),
    marital = factor(as.character(marital), levels = c("Married/CL", "Single", "Wid/Div/Sep")),
    race    = factor(as.character(race),    levels = c("White", "Black", "Puerto Rican", "Other")),
    site    = factor(as.character(site)),
    across(c(dm, chronic_htn, infant_sex), ~ factor(as.character(.x)))) %>%
  arrange(.imp, MOMID)
rm(w); invisible(gc())

imps <- split(d, d$.imp)
M <- length(imps)
d1 <- imps[[1]]
stopifnot(nrow(d1) == 40700, sum(d1$ih1) == 700, M == 50)

# ------------------------------------------------------------------------------
# 1. ARE MVM AND AI DIFFERENT VARIABLES?
# ------------------------------------------------------------------------------
cat("\nCross-tab of any MVM (rows) by any AI (columns), analytic sample:\n")
print(table(MVM = d1$mvm_any, AI = d1$ai_any))
cat("Agreement:", sprintf("%.1f%%", 100 * mean(d1$mvm_any == d1$ai_any)),
    " Correlation:", sprintf("%.3f", cor(d1$mvm_any, d1$ai_any)), "\n")
cat("IH by MVM: ", sum(d1$ih1[d1$mvm_any == 1]), "/", sum(d1$mvm_any == 1),
    " vs ", sum(d1$ih1[d1$mvm_any == 0]), "/", sum(d1$mvm_any == 0), "\n")
cat("IH by AI:  ", sum(d1$ih1[d1$ai_any == 1]), "/", sum(d1$ai_any == 1),
    " vs ", sum(d1$ih1[d1$ai_any == 0]), "/", sum(d1$ai_any == 0), "\n")

# ------------------------------------------------------------------------------
# 2-4. RE-FIT AND POOL
# ------------------------------------------------------------------------------
covs <- "age + bmi + race + educ + income + marital + smoking + parity + dm + chronic_htn + infant_sex"
models <- list(
  `Primary (covariates + site, weighted)` = list(rhs = paste(covs, "+ site"), wt = TRUE),
  `Secondary (covariates, weighted)`      = list(rhs = covs, wt = TRUE),
  `Covariates + site, unweighted`         = list(rhs = paste(covs, "+ site"), wt = FALSE),
  `Crude`                                 = list(rhs = NULL, wt = FALSE))

fit_one <- function(dd, expo, rhs, wt) {
  f <- as.formula(paste("ih1 ~", paste(c(expo, rhs), collapse = " + ")))
  dd$wt_ <- if (wt) dd$ipw else 1
  m <- geeglm(f, data = dd, family = poisson(link = "log"), id = MOMID,
              corstr = "exchangeable", weights = wt_)
  c(b = unname(coef(m)[expo]), v = unname(vcov(m)[expo, expo]))
}
rubin <- function(b, v) {
  qbar <- mean(b); W <- mean(v); B <- var(b); Tv <- W + (1 + 1 / length(b)) * B
  c(rr = exp(qbar), lci = exp(qbar - 1.96 * sqrt(Tv)), uci = exp(qbar + 1.96 * sqrt(Tv)))
}

res <- list(); by_imp <- list()
for (lab in names(models)) {
  for (expo in c("mvm_any", "ai_any")) {
    cat("Fitting", expo, "/", lab, "\n")
    fits <- t(sapply(imps, fit_one, expo = expo, rhs = models[[lab]]$rhs, wt = models[[lab]]$wt))
    p <- rubin(fits[, "b"], fits[, "v"])
    res[[length(res) + 1]] <- tibble(model = lab, exposure = expo,
      rr = p["rr"], lci = p["lci"], uci = p["uci"],
      rr_ci_4dp = sprintf("%.4f (%.4f, %.4f)", p["rr"], p["lci"], p["uci"]),
      rr_ci_2dp = sprintf("%.2f (%.2f, %.2f)", p["rr"], p["lci"], p["uci"]),
      cases = sum(d1$ih1), n = nrow(d1), n_imps = M)
    by_imp[[length(by_imp) + 1]] <- tibble(model = lab, exposure = expo, imp = seq_len(M),
                                           rr = exp(fits[, "b"]))
  }
}
res <- bind_rows(res); by_imp <- bind_rows(by_imp)

cat("\nPooled estimates (four decimals):\n")
print(res %>% select(model, exposure, rr_ci_4dp, rr_ci_2dp), n = Inf, width = Inf)

cat("\nPer-imputation primary aRR (min to max):\n")
print(by_imp %>% filter(grepl("^Primary", model)) %>% group_by(exposure) %>%
        summarise(min = round(min(rr), 4), max = round(max(rr), 4)))
cat("Imputations where the two primary aRRs are identical to four decimals:",
    by_imp %>% filter(grepl("^Primary", model)) %>%
      select(exposure, imp, rr) %>% pivot_wider(names_from = exposure, values_from = rr) %>%
      summarise(n = sum(round(mvm_any, 4) == round(ai_any, 4))) %>% pull(n), "of", M, "\n")

dir.create(out_dir, showWarnings = FALSE)
write_csv(res, file.path(out_dir, "verify_any_mvm_ai.csv"))
write_csv(by_imp, file.path(out_dir, "verify_any_mvm_ai_by_imp.csv"))
cat("Saved:", file.path(out_dir, "verify_any_mvm_ai.csv"), "\n")
