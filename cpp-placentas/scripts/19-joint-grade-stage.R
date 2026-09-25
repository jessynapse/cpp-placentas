# ==============================================================================
# 19-joint-grade-stage.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Dose-response with MVM and AI together, each in 3 categories:
#   MVM grade: none (0), low (score 1-2), high (score 3+)   [mvm3]
#   AI stage:  none, low, high (high = amnion or umbilical
#              artery involved, per Freedman et al. 2024)  [ai3]
#
#   1. 3 x 3 joint exposure (9 groups), reference = no MVM and no AI
#   2. Mutually adjusted: mvm3 + ai3 in the same model
#
# Same methods as 13-18: definite IH, 40,700 infants, 700 cases, 20
# imputations, Rubin's rules, modified Poisson GEE (exchangeable, MOMID).
# Models: crude, secondary (11 covariates + IPW), primary (+ site + IPW).
# Some joint cells are small. Estimates with fewer than 20 cases in a
# group are flagged.
#
# Outputs (outputs/):
#   joint_3x3_counts.csv   IH by MVM grade x AI stage
#   joint_3x3_models.csv   estimates for the 9 groups and the mutually
#                          adjusted model
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
counts <- d1 %>% group_by(joint, mvm3, ai3) %>%
  summarise(n_cases = sum(ih1 == 1), n_obs = n(), .groups = "drop") %>%
  mutate(ih_pct = fmt(n_cases, n_obs))
grid <- counts %>% select(mvm3, ai3, ih_pct) %>%
  pivot_wider(names_from = ai3, names_prefix = "AI ", values_from = ih_pct)

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

specs <- c("3x3 joint (ref: no MVM, no AI)" = "joint",
           "Mutually adjusted: MVM grade + AI stage" = "mvm3 + ai3")
res <- imap_dfr(specs, function(expo, sp) imap_dfr(models, function(extra, lab) {
  cat(" ", sp, "/", lab, "\n")
  run(expo, extra) %>% mutate(analysis = sp, model = lab)
})) %>%
  mutate(group = sub("^joint", "", term),
         group = recode(group, mvm3Low = "MVM low (adj. for AI stage)",
                        mvm3High = "MVM high (adj. for AI stage)",
                        ai3Low = "AI low (adj. for MVM grade)",
                        ai3High = "AI high (adj. for MVM grade)"),
         rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci)) %>%
  left_join(counts %>% transmute(group = as.character(joint), n_cases, n_obs, ih_pct),
            by = "group") %>%
  mutate(flag = ifelse(!is.na(n_cases) & n_cases < 20, "fewer than 20 cases", ""))

wide <- res %>% select(analysis, group, ih_pct, flag, model, rr_ci) %>%
  pivot_wider(names_from = model, values_from = rr_ci)

cat("\n==================== IH BY MVM GRADE (rows) x AI STAGE (columns) ====================\n")
print(grid, width = Inf)
cat("\n==================== ESTIMATES (definite IH) ====================\n")
print(wide, n = Inf, width = Inf)

write_csv(counts, file.path(out_dir, "joint_3x3_counts.csv"))
write_csv(wide,   file.path(out_dir, "joint_3x3_models.csv"))
cat("\nSaved: outputs/joint_3x3_counts.csv, outputs/joint_3x3_models.csv\n")
