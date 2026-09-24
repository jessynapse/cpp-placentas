# ==============================================================================
# 10-outcomemodels-corrected.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026. Corrected version of 05-outcomemodels.R
#
# Same models, covariates, GEE specification and Rubin's rules pooling as 05.
#
# CHANGES FROM 05:
#   1. Reads weighted_long_corrected.RDS from 09 (singletons only,
#      corrected chronic HTN, no covariates left unimputed)
#   2. plurality removed (all pregnancies are singletons)
#   3. Income reference = $4,000-$5,999, matching the manuscript
#   4. No p-values are computed or reported (CLAUDE.md rule)
#   5. Every model reports the number of IH cases and observations, and
#      categorical exposures report cases and observations per level
#      (outputs/corrected_group_counts.csv)
#   6. Site-adjusted models added for the V2 (suspect = present) outcome,
#      so every primary and V2 estimate has a with-site version
#   7. The 20 imputations are fitted in parallel (mclapply). Results are
#      identical to fitting them one at a time. On Windows set
#      n_cores <- 1.
#   8. Results saved to outputs/corrected_outcome_results.csv and to the
#      data folder as outcome_results_corrected.RDS
#   9. The Boston site_filter option in 05 is not carried over (unused)
#
# GEE Poisson (log link) with robust sandwich SE, exchangeable correlation,
# clustered on MOMID. IPW applied to the observed-outcome subset only.
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"   # folder with weighted_long_corrected.RDS
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
dir.create(out_dir, showWarnings = FALSE)
n_cores <- min(4, detectCores())

weighted <- readRDS(file.path(in_dir, "weighted_long_corrected.RDS"))

# ------------------------------------------------------------------------------
# FACTOR CODING (as 05, minus plurality, income reference changed)
# ------------------------------------------------------------------------------
weighted <- weighted %>%
  mutate(
    mvm3      = factor(mvm3,      levels = c("None","Low","High")),
    ai3       = factor(ai3,       levels = c("None","Low","High")),
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
  mutate(across(c(dm, chronic_htn, infant_sex, chorangioma,
                  mvm2, mvm_villous, mvm_vascular, ai), as.factor)) %>%
  mutate(
    ih_obs   = as.integer(!is.na(ih1) & exclusion == 0),
    mvm_num  = as.numeric(mvm),
    mvm2_num = as.numeric(as.character(mvm2)),
    ai_num   = as.numeric(as.character(ai))
  )

cat("Weighted long dims:", nrow(weighted), "x", ncol(weighted), "\n")
cat("Imputations:", max(weighted$.imp), "\n")
cat("Observed outcome rows per imp:",
    sum(weighted$ih_obs[weighted$.imp == 1] == 1), "\n")

# Split once so each parallel worker gets one imputation
imp_data <- weighted %>%
  filter(.imp > 0, ih_obs == 1) %>%
  arrange(MOMID) %>%
  split(.$.imp)

# ------------------------------------------------------------------------------
# RUBIN'S RULES POOLING FOR GEE (as 05, p-values removed)
# ------------------------------------------------------------------------------
pool_gee <- function(fits) {
  m        <- length(fits)
  coef_mat <- do.call(rbind, lapply(fits, coef))
  se_mat   <- do.call(rbind, lapply(fits, function(f) sqrt(diag(vcov(f)))))
  Q_bar <- colMeans(coef_mat)
  U_bar <- colMeans(se_mat^2)
  B     <- apply(coef_mat, 2, var)
  T_var <- U_bar + (1 + 1/m) * B
  SE    <- sqrt(T_var)
  tibble(
    term    = names(Q_bar),
    log_rr  = Q_bar,
    se      = SE,
    rr      = exp(Q_bar),
    lci     = exp(Q_bar - 1.96 * SE),
    uci     = exp(Q_bar + 1.96 * SE),
    n_imps  = m
  )
}

fmt_rr <- function(rr, lci, uci) sprintf("%.2f (%.2f, %.2f)", rr, lci, uci)

# ------------------------------------------------------------------------------
# COVARIATE SETS (as 05)
# ------------------------------------------------------------------------------
confounders_primary <- "age + bmi + race + educ + income + marital +
                        smoking + parity + dm + chronic_htn + infant_sex"
confounders_site    <- paste(confounders_primary, "+ site")

cov_vars_primary <- c("age", "bmi", "race", "educ", "income", "marital",
                      "smoking", "parity", "dm", "chronic_htn", "infant_sex")
cov_vars_site    <- c(cov_vars_primary, "site")

covariate_terms <- "^site|^race|^educ|^income|^marital|^smoking|^parity|^age$|^bmi$|^dm|^chronic_htn|^infant_sex"

# ------------------------------------------------------------------------------
# COUNTS: cases and observations, overall and per exposure level
# Covariates are complete after 08, so counts are identical across
# imputations. They are taken from imputation 1.
# ------------------------------------------------------------------------------
model_data <- function(d, covs, exp_vars, outcome) {
  d %>% filter(if_all(all_of(c(covs, exp_vars, outcome)), ~ !is.na(.)))
}

group_counts <- function(d, exp_vars, outcome, label, model) {
  out <- list()
  for (v in exp_vars) {
    x <- d[[v]]
    if (is.numeric(x) && length(unique(x)) > 3) next  # continuous
    out[[v]] <- d %>%
      group_by(level = as.character(.data[[v]])) %>%
      summarise(n_cases = sum(.data[[outcome]] == 1), n_obs = n(), .groups = "drop") %>%
      mutate(variable = v)
  }
  if (length(out) == 0) return(NULL)  # continuous exposure, no levels
  bind_rows(out) %>%
    mutate(exposure = label, outcome = outcome, model = model,
           ih_pct = sprintf("%.1f%% (%s/%s)", 100 * n_cases / n_obs,
                            format(n_cases, big.mark = ",", trim = TRUE),
                            format(n_obs, big.mark = ",", trim = TRUE)))
}

all_group_counts <- list()

# ------------------------------------------------------------------------------
# GENERIC MODEL RUNNER
# ------------------------------------------------------------------------------
run_gee <- function(exposure, label, outcome, model,
                    covs = character(0), conf_str = NULL, use_ipw = FALSE) {
  exp_vars <- trimws(unlist(strsplit(exposure, "\\*")))
  rhs <- if (is.null(conf_str)) exposure else paste(exposure, "+", conf_str)
  f   <- as.formula(paste(outcome, "~", rhs))

  fits <- mclapply(imp_data, function(d) {
    d <- model_data(d, covs, exp_vars, outcome) %>%
      mutate(across(where(is.factor), droplevels))
    # weights stored as a column: geeglm evaluates weights within data
    d$.w <- if (use_ipw) d$ipw else rep(1, nrow(d))
    tryCatch(geeglm(f, data = d, family = poisson(link = "log"),
                    id = MOMID, corstr = "exchangeable", weights = .w),
             error = function(e) e)
  }, mc.cores = n_cores)

  bad <- sapply(fits, function(x) !inherits(x, "geeglm"))
  if (any(bad)) message("  ", sum(bad), " fit(s) failed for ", label, " / ", model,
                        ": ", conditionMessage(fits[bad][[1]]))
  fits <- fits[!bad]
  if (length(fits) == 0) { warning("All fits failed: ", label); return(NULL) }

  d1 <- model_data(imp_data[[1]], covs, exp_vars, outcome)
  all_group_counts[[length(all_group_counts) + 1]] <<-
    group_counts(d1, exp_vars, outcome, label, model)

  pool_gee(fits) %>%
    filter(term != "(Intercept)", !grepl(covariate_terms, term)) %>%
    mutate(exposure = label, outcome = outcome, model = model,
           n_cases  = sum(d1[[outcome]] == 1),
           n_obs    = nrow(d1),
           rr_ci    = fmt_rr(rr, lci, uci))
}

# Wrappers
run_crude   <- function(e, l, o) run_gee(e, l, o, "crude")
run_adj     <- function(e, l, o, ipw)
  run_gee(e, l, o, ifelse(ipw, "adjusted+IPW", "adjusted"),
          covs = cov_vars_primary, conf_str = confounders_primary, use_ipw = ipw)
run_adj_site <- function(e, l, o, ipw)
  run_gee(e, l, o, ifelse(ipw, "adjusted+site+IPW", "adjusted+site"),
          covs = cov_vars_site, conf_str = confounders_site, use_ipw = ipw)

# ------------------------------------------------------------------------------
# EXPOSURE LISTS (as 05)
# ------------------------------------------------------------------------------
exposures_primary <- list(
  list(var = "mvm2", label = "MVM any vs none"),
  list(var = "mvm",  label = "MVM score (continuous)"),
  list(var = "ai",   label = "AI any vs none")
)
exposures_crude <- list(
  list(var = "mvm2_num", label = "MVM any vs none"),
  list(var = "mvm_num",  label = "MVM score (continuous)"),
  list(var = "ai_num",   label = "AI any vs none")
)
exposures_supp <- list(
  list(var = "mvm3",        label = "MVM grade (3-level)"),
  list(var = "ai3",         label = "AI stage (3-level)"),
  list(var = "mvm2 * ai",   label = "MVM x AI interaction"),
  list(var = "chorangioma", label = "Tumor (possible chorangioma)")
)

run_set <- function(exps, outcome, analysis, fun) {
  map_dfr(exps, function(e) { cat(" ", e$label, "\n"); fun(e$var, e$label, outcome) }) %>%
    mutate(analysis = analysis)
}

# ------------------------------------------------------------------------------
# RUN ALL MODELS
# ------------------------------------------------------------------------------
t0 <- Sys.time()
res <- list()
for (o in c("ih1", "ih2")) {
  ol <- ifelse(o == "ih1", "V1 definite", "V2 suspect = present")
  cat("\nCrude,", ol, "\n")
  res[[length(res)+1]] <- run_set(exposures_crude, o, paste(ol, "- crude"), run_crude)
  for (ipw in c(FALSE, TRUE)) {
    cat("\nPrimary exposures,", ol, ", no site, IPW =", ipw, "\n")
    res[[length(res)+1]] <- run_set(exposures_primary, o, paste(ol, "- primary"),
                                    function(e, l, oo) run_adj(e, l, oo, ipw))
    cat("\nPrimary exposures,", ol, ", with site, IPW =", ipw, "\n")
    res[[length(res)+1]] <- run_set(exposures_primary, o, paste(ol, "- site"),
                                    function(e, l, oo) run_adj_site(e, l, oo, ipw))
    cat("\nSupplemental exposures,", ol, ", no site, IPW =", ipw, "\n")
    res[[length(res)+1]] <- run_set(exposures_supp, o, paste(ol, "- supplemental"),
                                    function(e, l, oo) run_adj(e, l, oo, ipw))
  }
}
cat("\nModel fitting time:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

results <- bind_rows(res) %>%
  select(analysis, outcome, exposure, model, term, n_cases, n_obs,
         rr_ci, rr, lci, uci, log_rr, se, n_imps)
group_tab <- bind_rows(all_group_counts) %>%
  distinct() %>%
  select(exposure, outcome, model, variable, level, n_cases, n_obs, ih_pct)

# ------------------------------------------------------------------------------
# PRINT KEY RESULTS
# ------------------------------------------------------------------------------
show <- function(title, filt) {
  cat("\n======================================================================\n")
  cat(title, "\n")
  cat("======================================================================\n")
  results %>% filter({{ filt }}) %>%
    select(exposure, model, term, n_cases, n_obs, rr_ci) %>%
    print(n = Inf, width = Inf)
}
show("CRUDE, V1 definite",            analysis == "V1 definite - crude")
show("PRIMARY, V1, adjusted + IPW, no site vs with site",
     outcome == "ih1" & model %in% c("adjusted+IPW", "adjusted+site+IPW") &
       grepl("primary|site", analysis))
show("V2 (suspect = present), adjusted + IPW, no site vs with site",
     outcome == "ih2" & model %in% c("adjusted+IPW", "adjusted+site+IPW") &
       grepl("primary|site", analysis))
show("SUPPLEMENTAL, V1, adjusted + IPW",
     analysis == "V1 definite - supplemental" & model == "adjusted+IPW")

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
write_csv(results,   file.path(out_dir, "corrected_outcome_results.csv"))
write_csv(group_tab, file.path(out_dir, "corrected_group_counts.csv"))
saveRDS(results,     file.path(in_dir,  "outcome_results_corrected.RDS"))

cat("\nSaved: outputs/corrected_outcome_results.csv,",
    "outputs/corrected_group_counts.csv\n")
cat("Total result rows:", nrow(results), "\n")
