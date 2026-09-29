# ==============================================================================
# 31-maternal-fetal-ai.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Maternal vs fetal acute inflammation, requested by Linda Ernst ("I wonder
# if we should look at maternal vs fetal acute inflammation", "Also MIR and
# FIR?"). Supplemental table.
#
# Built from the neutrophil grades (0 to 3) in path2all.sas7bdat, the same
# 7 compartments that define AI_DI (see 14):
#   Maternal inflammatory response (MIR): amnion membrane roll (PA02_34),
#     chorion membrane roll (PA02_35), amnion placental surface (PA02_36),
#     chorion placental surface (PA02_37)
#   Fetal inflammatory response (FIR): umbilical vein (PA02_25), umbilical
#     artery (PA02_26), fetal surface vessels (PA02_38)
# A compartment is involved if its grade is above 0. Missing grades count
# as not involved (as in 14). Any MIR or FIR is identical to AI_DI (checked).
#
# Exposure (4 levels): no AI (reference), maternal only, fetal only,
# maternal and fetal.
#
# Corrected pipeline (07-10): 40,700 infants, 700 definite IH cases,
# 20 imputations, modified Poisson GEE (exchangeable, clustered on MOMID).
# Weights from 09 (denominator includes binary AI, not the MIR/FIR split).
# Models: crude, primary (11 covariates + site + IPW), secondary (11
# covariates + IPW, no site).
#
# Outputs (outputs/):
#   maternal_fetal_ai.csv          IH by group and RRs, one row per group
#   maternal_fetal_compartments.csv  compartments involved by group
# No p-values.
# ==============================================================================

library(tidyverse)
library(haven)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"   # folder with path2all.sas7bdat and weighted_long_corrected.RDS
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE)
n_cores <- min(4, detectCores())   # set to 1 on Windows

fmt <- function(x, n) sprintf("%.1f%% (%s/%s)", 100 * x / n,
                              format(x, big.mark = ",", trim = TRUE),
                              format(n, big.mark = ",", trim = TRUE))

# ------------------------------------------------------------------------------
# MATERNAL AND FETAL INFLAMMATION FROM COMPARTMENT GRADES
# ------------------------------------------------------------------------------
mat_comps <- c(amnion_roll = "PA02_34", chorion_roll = "PA02_35",
               amnion_surface = "PA02_36", chorion_surface = "PA02_37")
fet_comps <- c(umbilical_vein = "PA02_25", umbilical_artery = "PA02_26",
               fetal_surface_vessels = "PA02_38")
lv <- c("No AI", "Maternal only", "Fetal only", "Maternal and fetal")

path2 <- read_sas(file.path(in_dir, "path2all.sas7bdat"),
                  col_select = all_of(c("MOMID", "PREGID", "CHILDID",
                                        unname(c(mat_comps, fet_comps))))) %>%
  zap_labels() %>%
  mutate(across(all_of(unname(c(mat_comps, fet_comps))), ~ !is.na(.x) & .x > 0)) %>%
  rename(all_of(c(mat_comps, fet_comps))) %>%
  mutate(mir = if_any(all_of(names(mat_comps))),
         fir = if_any(all_of(names(fet_comps))),
         aimf = factor(case_when(!mir & !fir ~ lv[1], mir & !fir ~ lv[2],
                                 !mir & fir ~ lv[3], TRUE ~ lv[4]), levels = lv))

# ------------------------------------------------------------------------------
# DATA (as 16) + MIR/FIR
# ------------------------------------------------------------------------------
weighted <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS")) %>%
  left_join(path2, by = c("MOMID", "PREGID", "CHILDID"))

chk <- weighted %>% filter(.imp == 0)
stopifnot(!any(is.na(chk$aimf)))
stopifnot(all((chk$aimf != "No AI") == (as.numeric(as.character(chk$ai)) == 1)))
cat("Check passed: any MIR or FIR identical to AI_DI for all", nrow(chk),
    "pregnancies\n")

obs0 <- chk %>% filter(!is.na(ih1), exclusion == 0)
rm(chk)

# Which compartments are involved in each group (observed-outcome sample)
comps <- obs0 %>%
  filter(aimf != "No AI") %>%
  group_by(aimf) %>%
  summarise(n = n(), across(all_of(c(names(mat_comps), names(fet_comps))),
                            ~ sprintf("%.1f%%", 100 * mean(.x))),
            .groups = "drop")

imp_data <- weighted %>%
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
    site      = factor(site)
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), as.factor)) %>%
  select(.imp, MOMID, ih1, ipw, aimf, site, age, bmi, race, educ,
         income, marital, smoking, parity, dm, chronic_htn, infant_sex) %>%
  arrange(MOMID) %>%
  split(.$.imp)
rm(weighted); invisible(gc())
stopifnot(nrow(imp_data[[1]]) == 40700, sum(imp_data[[1]]$ih1) == 700)

# ------------------------------------------------------------------------------
# MODELS (as 16)
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
models <- c("Crude" = NA, "Primary: covariates + site + IPW" = "+ site",
            "Secondary: covariates + IPW, no site" = "")

run <- function(extra) {
  rhs <- if (is.na(extra)) "aimf" else paste("aimf +", confounders, extra)
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
  pool_gee(fits[!bad]) %>% filter(grepl("^aimf", term))
}

results <- imap_dfr(models, function(extra, lab) {
  cat(" ", lab, "\n")
  run(extra) %>% mutate(model = lab)
}) %>%
  mutate(group = sub("^aimf", "", term),
         rr_ci = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci))

counts <- imp_data[[1]] %>% group_by(aimf) %>%
  summarise(cases = sum(ih1 == 1), obs = n(), .groups = "drop") %>%
  transmute(group = as.character(aimf), `IH, % (cases/N)` = fmt(cases, obs))

table_out <- results %>%
  select(group, model, rr_ci) %>%
  pivot_wider(names_from = model, values_from = rr_ci) %>%
  bind_rows(tibble(group = "No AI", Crude = "1 REF",
                   `Primary: covariates + site + IPW` = "1 REF",
                   `Secondary: covariates + IPW, no site` = "1 REF"), .) %>%
  left_join(counts, by = "group") %>%
  select(`AI location` = group, `IH, % (cases/N)`, Crude,
         Primary = `Primary: covariates + site + IPW`,
         Secondary = `Secondary: covariates + IPW, no site`) %>%
  mutate(n_cases = 700, n_obs = 40700)

cat("\n==================== COMPARTMENTS INVOLVED BY GROUP ====================\n")
print(comps, width = Inf)
cat("\n==================== MATERNAL VS FETAL AI (definite IH) ====================\n")
print(table_out, width = Inf)

write_csv(table_out, file.path(out_dir, "maternal_fetal_ai.csv"))
write_csv(comps, file.path(out_dir, "maternal_fetal_compartments.csv"))
cat("\nSaved: outputs/maternal_fetal_ai.csv, outputs/maternal_fetal_compartments.csv\n")
