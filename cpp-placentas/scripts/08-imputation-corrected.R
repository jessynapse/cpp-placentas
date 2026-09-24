# ==============================================================================
# 08-imputation-corrected.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026. Corrected copy of 03-imputation.R
#
# CHANGES FROM 03 (each marked "# CHANGED:" below):
#   1. Reads analytic_sample_corrected.RDS from 07 (singletons only,
#      corrected chronic HTN)
#   2. gest_age and birthweight removed from the imputation model. 03
#      treated them as fully observed predictors, but 225 women are
#      missing gest_age and 56 are missing birthweight. mice cannot
#      impute a row whose predictor is missing and not imputed, so 57
#      women kept missing covariates in every imputation and were
#      silently dropped from the adjusted models in 05. Neither variable
#      is a model covariate (both are potential mediators), so they are
#      no longer used as predictors. They are added back to the long
#      dataset afterwards, unimputed, like the other perinatal variables.
#   3. plurality removed from the imputation model. After 07 every
#      pregnancy is a singleton, so it is constant.
#   4. The 20 imputations run in 4 parallel chains of 5 (seeds 12345 to
#      12348), combined with ibind(). Same model, about 4x faster.
#      Results differ from a single chain only through the random draws.
#      mclapply forks on macOS and Linux. On Windows set n_cores <- 1.
#   5. ai and ai3 are temporarily imputed. 03 treated them as fully
#      observed, but 2 placentas are missing AI, and one of those women
#      is also missing income, education and infant sex, so mice could
#      not impute her. After imputation, ai and ai3 are reset to their
#      observed values (NA for those 2), so the exposure is never
#      imputed and AI models still exclude them.
#   6. A check stops the script if any covariate is still missing after
#      imputation.
#   7. Outputs saved with a _corrected suffix.
#
# Original 03 notes:
# Imputes missing covariates only — outcomes (ih1, ih2) never imputed
# Exposures (mvm, ai) fully observed — included as predictors only
# m = 20 imputations, maxit = 10, seed = 12345
#
# Parity missingness corrected via na_tag() in 01-load-clean.R
# Expected parity missingness now ~1% (down from 28.6%)
# Logged events: MICE drops redundant exposure predictors
# (mvm2 collinear with mvm3, ai collinear with ai3) — imputation valid
# Density plot skipped — dm/chronic_htn too sparse for kernel density
# ==============================================================================

library(tidyverse)
library(mice)
library(VIM)

setwd("/Users/wongjj/Downloads")

# CHANGED: corrected analytic sample from 07
df <- readRDS("analytic_sample_corrected.RDS")

# ------------------------------------------------------------------------------
# FACTOR CODING
# ------------------------------------------------------------------------------
df2 <- df %>%
  mutate(
    mvm3        = factor(mvm3,        levels = 0:2,
                         labels = c("None","Low","High")),
    ai3         = factor(ai3,         levels = 0:2,
                         labels = c("None","Low","High")),
    smoking     = factor(smoking,     levels = 0:2,
                         labels = c("Non-smoker","1-19/day","20+/day")),
    parity      = factor(parity,      levels = 0:4,
                         labels = c("0","1","2","3","4+")),
    plurality   = factor(plurality,   levels = 1:2,
                         labels = c("Singleton","Multiple")),
    income      = factor(income,      levels = 1:7,
                         labels = c("<=$1,999","$2,000-$3,999",
                                    "$4,000-$5,999","$6,000-$7,999",
                                    "$8,000-$9,999","$10,000-$14,999",
                                    ">=$15,000")),
    educ        = factor(educ,        levels = 0:3,
                         labels = c("<HS","Some HS","HS grad","College+")),
    marital     = factor(marital,     levels = 0:2,
                         labels = c("Married/CL","Single","Wid/Div/Sep")),
    race        = factor(race,        levels = c(1,2,4,8),
                         labels = c("White","Black","Puerto Rican","Other")),
    infant_sex  = factor(infant_sex,  levels = 0:1,
                         labels = c("Female","Male")),
    dm          = factor(dm,          levels = 0:1, labels = c("No","Yes")),
    chronic_htn = factor(chronic_htn, levels = 0:1, labels = c("No","Yes")),
    site        = factor(site)
  )

# ------------------------------------------------------------------------------
# MISSINGNESS CHECK
# ------------------------------------------------------------------------------
cat("Analytic sample N:", nrow(df2), "\n")
cat("\nMissingness (covariates only):\n")
# CHANGED: plurality dropped (constant after 07)
miss_vars <- c("age","bmi","smoking","parity","income",
               "educ","marital","infant_sex","race","dm","chronic_htn","site")
for (v in miss_vars) {
  n <- sum(is.na(df2[[v]]))
  cat(sprintf("  %-15s %5d (%s%%)\n", v, n, round(100*n/nrow(df2), 1)))
}

cat("\nParity check — should be ~1% missing after na_tag() fix:\n")
cat("  Parity missing:  ", sum(is.na(df2$parity)), "\n")
cat("  Parity pct miss: ", round(100*mean(is.na(df2$parity)), 1), "%\n")
cat("  Parity = 0 (primigravida): ", sum(df2$parity == "0", na.rm=TRUE), "\n")

pdf("mice_missing_corrected.pdf", width = 12, height = 6)
aggr(df2 %>% select(all_of(miss_vars)),
     col = c("navyblue","red"), numbers = TRUE,
     sortVars = TRUE, cex.axis = 0.6,
     ylab = c("Missing data","Pattern"))
dev.off()

# ------------------------------------------------------------------------------
# PREPARE MICE INPUT
# Excluded: IDs, ih1, ih2, exclusion, perinatal outcomes
# Included: all covariates + exposures as predictors
# ------------------------------------------------------------------------------
vars_to_impute <- c(
  "mvm", "mvm3", "mvm2", "ai", "ai3",      # exposures — predictors only
  "age", "bmi", "smoking", "parity",
  "income", "educ", "marital",
  "infant_sex", "race", "dm", "chronic_htn",
  "site")  # CHANGED: plurality, gest_age, birthweight removed

df_mi <- df2 %>% select(all_of(vars_to_impute))

cat("\nVariables entering MICE:", ncol(df_mi), "\n")
cat("Rows:", nrow(df_mi), "\n")

# ------------------------------------------------------------------------------
# INITIALIZE + SET METHODS
# ------------------------------------------------------------------------------
ini  <- mice(df_mi, maxit = 0, seed = 12345)
meth <- ini$method
pred <- ini$predictorMatrix

# Fully observed — set to "" so MICE does not impute them
# CHANGED: gest_age, birthweight, plurality no longer in the model
# CHANGED: ai and ai3 removed from fully_obs (2 missing). They are
# imputed only so other covariates can be imputed, then reset below.
fully_obs <- c("mvm","mvm3","mvm2","age","race","site")
meth[fully_obs] <- ""
pred[fully_obs, ] <- 0

cat("\nMICE imputation methods:\n")
print(meth[meth != ""])

# ------------------------------------------------------------------------------
# RUN MICE
# m = 20, maxit = 10
# CHANGED: 4 parallel chains of m = 5 (seeds 12345 to 12348), then ibind()
# ------------------------------------------------------------------------------
library(parallel)
n_cores <- min(4, detectCores())
cat("\nRunning MICE: m=20 (4 chains x 5), maxit=10, cores =", n_cores, "...\n")

imp_list <- mclapply(1:4, function(k) {
  mice(df_mi,
       m               = 5,
       method          = meth,
       predictorMatrix = pred,
       seed            = 12344 + k,
       maxit           = 10,
       printFlag       = FALSE)
}, mc.cores = n_cores)

failed <- sapply(imp_list, function(x) !inherits(x, "mids"))
if (any(failed)) stop("MICE chain(s) failed: ", paste(which(failed), collapse = ", "),
                      "\n", paste(sapply(imp_list[failed], as.character), collapse = "\n"))

df_imp <- Reduce(ibind, imp_list)
cat("Combined imputations:", df_imp$m, "\n")

# ------------------------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS
# Density plot skipped — dm/chronic_htn too sparse for kernel density
# ------------------------------------------------------------------------------
cat("\nSaving convergence plot...\n")
pdf("mice_convergence_corrected.pdf", width = 12, height = 8)
plot(df_imp)
dev.off()

pdf("mice_density_corrected.pdf", width = 12, height = 8)
tryCatch(densityplot(df_imp), error = function(e) {
  message("Density plot skipped: ", e$message)
  stripplot(df_imp, pch = 20, cex = 0.8)
})
dev.off()

# ------------------------------------------------------------------------------
# BUILD LONG DATASET
# Add back IDs, outcomes, and perinatal variables (never imputed)
# ih1/ih2 missingness preserved — handled downstream by IPW
# ------------------------------------------------------------------------------
long <- complete(df_imp, "long", include = TRUE)

long <- long %>%
  mutate(
    MOMID        = rep(df2$MOMID,        times = df_imp$m + 1),
    PREGID       = rep(df2$PREGID,       times = df_imp$m + 1),
    CHILDID      = rep(df2$CHILDID,      times = df_imp$m + 1),
    ih1          = rep(df2$ih1,          times = df_imp$m + 1),
    ih2          = rep(df2$ih2,          times = df_imp$m + 1),
    exclusion    = rep(df2$exclusion,    times = df_imp$m + 1),
    preterm      = rep(df2$preterm,      times = df_imp$m + 1),
    preterm3     = rep(df2$preterm3,     times = df_imp$m + 1),
    SGA          = rep(df2$SGA,          times = df_imp$m + 1),
    bw_for_ga_z  = rep(df2$bw_for_ga_z, times = df_imp$m + 1),
    labor_onset  = rep(df2$labor_onset,  times = df_imp$m + 1),
    chorangioma  = rep(df2$chorangioma,  times = df_imp$m + 1),
    mvm_villous  = rep(df2$mvm_villous,  times = df_imp$m + 1),
    mvm_vascular = rep(df2$mvm_vascular, times = df_imp$m + 1),
    # CHANGED: exposure reset to observed values, never imputed
    ai           = rep(df2$ai,           times = df_imp$m + 1),
    ai3          = rep(df2$ai3,          times = df_imp$m + 1),
    # CHANGED: added back unimputed, no longer in the imputation model
    gest_age     = rep(df2$gest_age,     times = df_imp$m + 1),
    birthweight  = rep(df2$birthweight,  times = df_imp$m + 1)
  )

# ------------------------------------------------------------------------------
# VERIFY
# ------------------------------------------------------------------------------
cat("\nVerification:\n")
cat("  Long dims:          ", nrow(long), "x", ncol(long), "\n")
cat("  Imputations:        ", max(long$.imp), "\n")
cat("  IH cases imp 0:     ", sum(long$ih1[long$.imp==0]==1, na.rm=TRUE), "\n")
cat("  IH cases imp 1:     ", sum(long$ih1[long$.imp==1]==1, na.rm=TRUE), "\n")
cat("  Parity miss imp 0:  ", sum(is.na(long$parity[long$.imp==0])), "\n")
cat("  Parity miss imp 1:  ", sum(is.na(long$parity[long$.imp==1])), "\n")

# CHANGED: stop if any covariate is still missing after imputation
cov_check <- c("age","bmi","smoking","parity","income","educ","marital",
               "infant_sex","race","dm","chronic_htn","site")
n_left <- sapply(cov_check, function(v) sum(is.na(long[[v]][long$.imp > 0])))
cat("  Covariate values still missing across imputations 1-20:", sum(n_left), "\n")
if (sum(n_left) > 0) {
  print(n_left[n_left > 0])
  stop("Imputation left missing covariate values")
}

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
# CHANGED: _corrected file names
saveRDS(df_imp, "mice_obj_corrected.RDS")
saveRDS(long,   "imputed_long_corrected.RDS")
write_csv(long, "imputed_long_corrected.csv", na = "")

cat("\nSaved: mice_obj_corrected.RDS, imputed_long_corrected.RDS, imputed_long_corrected.csv\n")
cat("mice_convergence_corrected.pdf, mice_density_corrected.pdf\n")
cat("N:", nrow(df2), "| Imputations: 20 | Rows in long:", nrow(long), "\n")