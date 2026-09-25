# ==============================================================================
# 15-pathologist.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Does adjustment for the individual placental examiner (pathologist)
# explain the attenuation seen with study site adjustment?
#
# Motivation: Freedman et al. 2024 (ATVB) restricted CPP sibling analyses to
# placentas read by the same pathologist because of interobserver
# variability. path1all.sas7bdat has an examiner code (PA101_15).
#
# Examiner definition: examiner codes repeat across sites (37 of 66 codes
# appear at more than one site), so they are local numbers within a site.
# A pathologist is therefore site x PA101_15. Examiners with < 100
# placentas in the analytic sample, and missing codes (n = 74), are
# grouped into one "other" category per site. This gives 62 named
# examiners (90.7% of placentas) plus up to 12 site "other" groups.
# Because examiners are nested within site, adjusting for examiner also
# adjusts for site.
#
# CAUTION: PA101_15 is on the gross examination form (PA101). It is not
# certain that the same person read the histology (PA02 form). Confirm
# with Alexa or Linda before describing it as the histology pathologist.
#
# Models (definite IH, 40,700 infants, 700 cases, 20 imputations, IPW
# from 09, Rubin's rules), for MVM any, MVM score, AI any:
#   Secondary:  11 covariates + IPW (no site)
#   Primary:    11 covariates + site + IPW
#   Examiner:   11 covariates + examiner (site x code) + IPW
#
# Outputs (outputs/):
#   pathologist_models.csv        estimates
#   pathologist_by_site.csv       within-site spread of MVM and AI by examiner
# No p-values.
# ==============================================================================

library(tidyverse)
library(haven)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"   # folder with path1all.sas7bdat and weighted_long_corrected.RDS
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE)
n_cores <- min(4, detectCores())   # set to 1 on Windows

# ------------------------------------------------------------------------------
# EXAMINER VARIABLE
# ------------------------------------------------------------------------------
p1 <- read_sas(file.path(in_dir, "path1all.sas7bdat"),
               col_select = c(MOMID, PREGID, CHILDID, PA101_15)) %>%
  zap_labels()

weighted <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  left_join(p1, by = c("MOMID", "PREGID", "CHILDID"))

base <- weighted %>% filter(.imp == 0)
ex_n <- base %>% count(site, PA101_15, name = "n_placentas")
weighted <- weighted %>%
  left_join(ex_n, by = c("site", "PA101_15")) %>%
  mutate(examiner = ifelse(!is.na(PA101_15) & n_placentas >= 100,
                           paste0("site", site, "_ex", PA101_15),
                           paste0("site", site, "_other")))

# Within-site spread (descriptive)
by_site <- weighted %>% filter(.imp == 0) %>%
  group_by(site, examiner) %>%
  summarise(n = n(), mvm_pct = 100 * mean(mvm2 == 1), ai_pct = 100 * mean(ai == 1),
            .groups = "drop") %>%
  filter(!grepl("other", examiner)) %>%
  group_by(site) %>%
  summarise(examiners_100plus = n(), placentas = sum(n),
            any_mvm_range = sprintf("%.0f%% to %.0f%%", min(mvm_pct), max(mvm_pct)),
            any_ai_range  = sprintf("%.0f%% to %.0f%%", min(ai_pct),  max(ai_pct)),
            .groups = "drop")
cat("Examiner groups:", n_distinct(weighted$examiner), "\n")
print(by_site, n = Inf)

# ------------------------------------------------------------------------------
# DATA CODING (as 10)
# ------------------------------------------------------------------------------
weighted <- weighted %>%
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
    examiner  = factor(examiner),
    ih_obs    = as.integer(!is.na(ih1) & exclusion == 0),
    mvm2      = factor(mvm2),
    ai        = factor(ai)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor))

imp_data <- weighted %>% filter(.imp > 0, ih_obs == 1) %>% arrange(MOMID) %>% split(.$.imp)

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
         lci = exp(Q_bar - 1.96 * SE), uci = exp(Q_bar + 1.96 * SE), n_imps = m)
}

confounders <- "age + bmi + race + educ + income + marital +
                smoking + parity + dm + chronic_htn + infant_sex"
models <- c("Secondary: covariates, no site" = "",
            "Primary: covariates + site"     = "+ site",
            "Covariates + examiner"          = "+ examiner")
exposures <- c(mvm2 = "MVM any vs none", mvm = "MVM score, per 1 point",
               ai = "AI any vs none")
terms <- c(mvm2 = "mvm21", mvm = "mvm", ai = "ai1")

run <- function(e, extra) {
  f <- as.formula(paste("ih1 ~", e, "+", confounders, extra))
  fit_one <- function(d) {
    d <- d %>% mutate(examiner = droplevels(examiner))
    tryCatch(geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                    corstr = "exchangeable", weights = ipw),
             error = function(err) err)
  }
  # the examiner model is memory heavy, so it uses fewer parallel workers
  cores <- if (grepl("examiner", extra)) 2 else n_cores
  fits <- mclapply(imp_data, fit_one, mc.cores = cores)
  # a crashed worker returns NULL or try-error: refit those one at a time
  redo <- which(!sapply(fits, function(x) inherits(x, "geeglm")))
  for (i in redo) fits[[i]] <- fit_one(imp_data[[i]])
  bad <- sapply(fits, function(x) !inherits(x, "geeglm"))
  if (any(bad)) message("  ", sum(bad), " fit(s) failed after retry: ", e, " ", extra)
  pool_gee(fits[!bad]) %>% filter(term == terms[[e]])
}

results <- map_dfr(names(exposures), function(e) imap_dfr(models, function(extra, lab) {
  cat(" ", exposures[[e]], "/", lab, "\n")
  run(e, extra) %>% mutate(exposure = exposures[[e]], model = lab)
})) %>%
  mutate(n_cases = sum(imp_data[[1]]$ih1 == 1), n_obs = nrow(imp_data[[1]]),
         rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci)) %>%
  select(exposure, model, n_cases, n_obs, rr_ci, rr, lci, uci, n_imps)

cat("\n==================== EXAMINER ADJUSTMENT (adjusted + IPW) ====================\n")
print(results %>% select(exposure, model, n_cases, n_obs, rr_ci), n = Inf, width = Inf)

write_csv(results, file.path(out_dir, "pathologist_models.csv"))
write_csv(by_site, file.path(out_dir, "pathologist_by_site.csv"))
cat("\nSaved: outputs/pathologist_models.csv, pathologist_by_site.csv\n")
