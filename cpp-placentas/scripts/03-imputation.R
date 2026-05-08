# ==============================================================================
# 03-imputation.R
# CPP Placental Pathology → Infantile Hemangioma
# Updated May 2026
#
# Imputes missing covariates only — outcomes (ih1, ih2) never imputed
# Exposures (mvm, ai) fully observed — included as predictors only
# m = 20 imputations, maxit = 10, seed = 12345
#
# Logged events (1800): MICE dropped mvm2/ai as redundant predictors
# due to collinearity with mvm3/ai3 — expected, imputation valid
# Density plot skipped — dm/chronic_htn too sparse for kernel density
# ==============================================================================

library(tidyverse)
library(mice)
library(VIM)

setwd("/Users/wongjj/Downloads")

df <- readRDS("analytic_sample.RDS")

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
miss_vars <- c("age","bmi","smoking","parity","plurality","income",
               "educ","marital","infant_sex","race","dm","chronic_htn","site")
for (v in miss_vars) {
  n <- sum(is.na(df2[[v]]))
  cat(sprintf("  %-15s %5d (%s%%)\n", v, n, round(100*n/nrow(df2), 1)))
}

pdf("mice_missing.pdf", width = 12, height = 6)
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
  "plurality", "income", "educ", "marital",
  "infant_sex", "race", "dm", "chronic_htn",
  "site", "gest_age", "birthweight"
)

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
fully_obs <- c("mvm","mvm3","mvm2","ai","ai3","age",
               "gest_age","birthweight","race","site","plurality")
meth[fully_obs] <- ""
pred[fully_obs, ] <- 0

cat("\nMICE imputation methods:\n")
print(meth[meth != ""])

# ------------------------------------------------------------------------------
# RUN MICE
# m = 20, maxit = 10, seed = 12345
# Note: 1800 logged events expected — MICE drops redundant exposure predictors
# (mvm2 collinear with mvm3, ai collinear with ai3) — imputation valid
# ------------------------------------------------------------------------------
cat("\nRunning MICE: m=20, maxit=10...\n")

df_imp <- mice(df_mi,
               m               = 20,
               method          = meth,
               predictorMatrix = pred,
               seed            = 12345,
               maxit           = 10,
               printFlag       = TRUE)

# ------------------------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS
# Density plot skipped — dm/chronic_htn too sparse for kernel density
# ------------------------------------------------------------------------------
cat("\nSaving convergence plot...\n")
pdf("mice_convergence.pdf", width = 12, height = 8)
plot(df_imp)
dev.off()

pdf("mice_density.pdf", width = 12, height = 8)
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
    MOMID        = rep(df2$MOMID,        times = 21),
    PREGID       = rep(df2$PREGID,       times = 21),
    CHILDID      = rep(df2$CHILDID,      times = 21),
    ih1          = rep(df2$ih1,          times = 21),
    ih2          = rep(df2$ih2,          times = 21),
    exclusion    = rep(df2$exclusion,    times = 21),
    preterm      = rep(df2$preterm,      times = 21),
    preterm3     = rep(df2$preterm3,     times = 21),
    SGA          = rep(df2$SGA,          times = 21),
    bw_for_ga_z  = rep(df2$bw_for_ga_z, times = 21),
    labor_onset  = rep(df2$labor_onset,  times = 21),
    chorangioma  = rep(df2$chorangioma,  times = 21),
    mvm_villous  = rep(df2$mvm_villous,  times = 21),
    mvm_vascular = rep(df2$mvm_vascular, times = 21)
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

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
saveRDS(df_imp, "mice_obj.RDS")
saveRDS(long,   "imputed_long.RDS")
write_csv(long, "imputed_long.csv", na = "")

cat("\nSaved: mice_obj.RDS, imputed_long.RDS, imputed_long.csv\n")
cat("mice_convergence.pdf, mice_density.pdf\n")
cat("N:", nrow(df2), "| Imputations: 20 | Rows in long:", nrow(long), "\n")