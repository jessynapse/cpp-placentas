# ==============================================================================
# 37-models-m50.R
# CPP Placental Pathology -> Infantile Hemangioma
# October 2026
#
# Refits every model reported in the manuscript on the 50 imputed datasets
# from 35 and 36, as required by Paediatric and Perinatal Epidemiology
# (minimum of 50 imputations), and adds E-values (VanderWeele and Ding
# 2017), also required by PPE.
#
# Same specifications as the 20-imputation scripts it replaces:
#   10  suspect-IH sensitivity (ih2)
#   13, 17  sequential adjustment and Boston stratification
#   14  four exposures (any MVM, MVM score, any AI, AI score)
#   16, 18  joint exposure and mutual adjustment
#   20  MVM grade and AI stage
#   21  site as a random effect
#   22  site terms from the primary model
#   29  main Table 2
#   30  RERI, attributable proportion, synergy index, ratio of RRs
#   31  AI by maternal or foetal involvement
#
# Modified Poisson GEE (exchangeable, clustered on MOMID), Rubin's rules.
# Models: crude, primary (11 covariates + site + IPW), secondary
# (11 covariates + IPW). Definite IH unless stated.
#
# Outputs (outputs/, all with an _m50 suffix and the same layout as the
# 20-imputation files, so the manuscript builder can switch to them):
#   table2_main_m50.csv, synergy_reri_m50.csv, mutual_4exp_m50.csv,
#   grade_stage_separate_m50.csv, maternal_fetal_ai_m50.csv,
#   seq_cumulative_4exp_m50.csv, boston_4exp_m50.csv,
#   site_random_effect_models_m50.csv, ih2_m50.csv,
#   site_table_percentages_m50.csv, evalues_m50.csv, all_estimates_m50.csv
# No p-values.
# ==============================================================================

library(tidyverse)
library(haven)
library(geepack)
library(lme4)
library(parallel)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
n_cores <- min(4, detectCores())   # set to 1 on Windows

fmt_ci <- function(rr, l, u) sprintf("%.2f (%.2f, %.2f)", rr, l, u)
fmt_n  <- function(x) format(x, big.mark = ",", trim = TRUE)
pct    <- function(k, n) sprintf("%.1f%% (%s/%s)", 100 * k / n, fmt_n(k), fmt_n(n))

# ------------------------------------------------------------------------------
# DATA
# ------------------------------------------------------------------------------
mat_c <- c("PA02_34", "PA02_35", "PA02_36", "PA02_37")
fet_c <- c("PA02_25", "PA02_26", "PA02_38")
path2 <- read_sas(file.path(in_dir, "path2all.sas7bdat"),
                  col_select = all_of(c("MOMID", "PREGID", "CHILDID", mat_c, fet_c))) %>%
  zap_labels() %>%
  mutate(across(all_of(c(mat_c, fet_c)), ~ !is.na(.x) & .x > 0),
         ai_count = rowSums(across(all_of(c(mat_c, fet_c)))),
         mir = if_any(all_of(mat_c)), fir = if_any(all_of(fet_c)),
         aimf = factor(case_when(!mir & !fir ~ "No AI", mir & !fir ~ "Maternal only",
                                 !mir & fir ~ "Fetal only", TRUE ~ "Maternal and fetal"),
                       levels = c("No AI", "Maternal only", "Fetal only", "Maternal and fetal"))) %>%
  select(MOMID, PREGID, CHILDID, ai_count, aimf)

# For a quick code test only: M50_TEST=2 uses the 20-imputation file and
# keeps 2 imputations, writing outputs with a _test suffix.
test_n <- as.integer(Sys.getenv("M50_TEST", "0"))
w <- readRDS(file.path(in_dir, if (test_n > 0) "weighted_long_corrected.RDS" else "weighted_long_m50.RDS"))
if (test_n > 0) w <- w[w$.imp <= test_n, ]
M <- max(w$.imp)
stopifnot(test_n > 0 || M >= 50)
sfx <- if (test_n > 0) "_test" else "_m50"

imp <- w %>%
  filter(.imp > 0, !is.na(ih1), exclusion == 0) %>%
  left_join(path2, by = c("MOMID", "PREGID", "CHILDID")) %>%
  mutate(
    smoking = factor(as.character(smoking), levels = c("Non-smoker", "1-19/day", "20+/day")),
    parity  = factor(as.character(parity),  levels = c("0", "1", "2", "3", "4+")),
    income  = factor(as.character(income),
                     levels = c("$4,000-$5,999", "<=$1,999", "$2,000-$3,999", "$6,000-$7,999",
                                "$8,000-$9,999", "$10,000-$14,999", ">=$15,000")),
    educ    = factor(as.character(educ),    levels = c("<HS", "Some HS", "HS grad", "College+")),
    marital = factor(as.character(marital), levels = c("Married/CL", "Single", "Wid/Div/Sep")),
    race    = factor(as.character(race),    levels = c("White", "Black", "Puerto Rican", "Other")),
    site    = factor(as.character(site),
                     levels = c("5", "10", "15", "31", "37", "45", "50", "55", "60", "66", "71", "82")),
    mvm2n   = as.numeric(as.character(mvm2)),
    ain     = as.numeric(as.character(ai)),
    mvm     = as.numeric(mvm),
    mvm3    = factor(as.character(mvm3), levels = c("None", "Low", "High")),
    ai3     = factor(as.character(ai3),  levels = c("None", "Low", "High")),
    joint   = factor(case_when(mvm2n == 0 & ain == 0 ~ "Neither", mvm2n == 1 & ain == 0 ~ "MVM only",
                               mvm2n == 0 & ain == 1 ~ "AI only", TRUE ~ "Both MVM and AI"),
                     levels = c("Neither", "MVM only", "AI only", "Both MVM and AI")),
    age_c   = as.numeric(scale(age)), bmi_c = as.numeric(scale(bmi))
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), ~ factor(as.character(.x)))) %>%
  select(.imp, MOMID, ih1, ih2, ipw, site, age, bmi, age_c, bmi_c, race, educ, income,
         marital, smoking, parity, dm, chronic_htn, infant_sex, mvm2n, mvm, ain, ai_count,
         mvm3, ai3, joint, aimf) %>%
  arrange(MOMID) %>%
  split(.$.imp)
rm(w); invisible(gc())
d1 <- imp[[1]]
stopifnot(nrow(d1) == 40700, sum(d1$ih1) == 700, length(imp) == M)
cat("Imputations:", length(imp), "\n")

# ------------------------------------------------------------------------------
# FIT AND POOL (Rubin's rules, with pooled covariance for the RERI)
# ------------------------------------------------------------------------------
covs <- c("age", "bmi", "race", "educ", "income", "marital", "smoking", "parity",
          "dm", "chronic_htn", "infant_sex")
spec <- list(Crude = list(covs = character(0), w = FALSE),
             Primary = list(covs = c(covs, "site"), w = TRUE),
             Secondary = list(covs = covs, w = TRUE))

fit_pool <- function(expo, covs, use_w, outcome = "ih1", data = imp, keep = NULL) {
  rhs <- paste(c(expo, covs), collapse = " + ")
  f <- as.formula(paste(outcome, "~", rhs))
  one <- function(d) {
    d$.w <- if (use_w) d$ipw else rep(1, nrow(d))
    tryCatch({
      m <- geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                  corstr = "exchangeable", weights = .w)
      list(b = coef(m), V = vcov(m))
    }, error = function(e) NULL)
  }
  fits <- mclapply(data, one, mc.cores = n_cores)
  for (i in which(sapply(fits, is.null))) fits[[i]] <- one(data[[i]])
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) < length(data)) message("  ", length(data) - length(fits), " fit(s) failed: ", rhs)
  m <- length(fits)
  B <- do.call(rbind, lapply(fits, `[[`, "b"))
  terms <- colnames(B)
  if (!is.null(keep)) terms <- grep(keep, terms, value = TRUE)
  B <- B[, terms, drop = FALSE]
  W <- Reduce(`+`, lapply(fits, function(x) x$V[terms, terms, drop = FALSE])) / m
  V <- W + (1 + 1/m) * cov(B)
  b <- colMeans(B); se <- sqrt(diag(V))
  list(est = tibble(term = terms, b = b, se = se, rr = exp(b),
                    lci = exp(b - 1.96 * se), uci = exp(b + 1.96 * se),
                    rr_ci = fmt_ci(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)), n_imps = m),
       V = V)
}

run3 <- function(expo, keep, outcome = "ih1", data = imp) {
  imap(spec, function(s, lab) { cat("   ", expo, "/", lab, "\n"); fit_pool(expo, s$covs, s$w, outcome, data, keep) })
}

all_est <- list()
add_est <- function(res, analysis) {
  all_est[[length(all_est) + 1]] <<- imap_dfr(res, ~ .x$est %>% mutate(model = .y)) %>% mutate(analysis = analysis)
}

# ------------------------------------------------------------------------------
# 1. FOUR EXPOSURES, JOINT, MUTUAL, GRADE/STAGE, MATERNAL/FOETAL
# ------------------------------------------------------------------------------
cat("\n== Four exposures\n")
r_mvm2 <- run3("mvm2n", "^mvm2n$|^site");  add_est(r_mvm2, "Any MVM")
r_mvm  <- run3("mvm", "^mvm$");            add_est(r_mvm, "MVM score")
r_ai   <- run3("ain", "^ain$|^site");      add_est(r_ai, "Any AI")
r_aic  <- run3("ai_count", "^ai_count$");  add_est(r_aic, "AI score")
cat("\n== Joint\n");  r_joint <- run3("joint", "^joint"); add_est(r_joint, "Joint")
cat("\n== Mutual\n")
r_mut  <- run3("mvm2n + ain", "^mvm2n$|^ain$");       add_est(r_mut, "Mutual binary")
r_mutc <- run3("mvm + ai_count", "^mvm$|^ai_count$"); add_est(r_mutc, "Mutual score")
cat("\n== Grade and stage\n")
r_g <- run3("mvm3", "^mvm3"); add_est(r_g, "MVM grade")
r_s <- run3("ai3", "^ai3");   add_est(r_s, "AI stage")
cat("\n== Maternal or foetal\n"); r_mf <- run3("aimf", "^aimf"); add_est(r_mf, "AI involvement")

get <- function(res, model, term) res[[model]]$est %>% filter(.data$term == !!term) %>% pull(rr_ci)

# Table 2
ih <- function(x, lvl) { k <- d1$ih1[x == lvl]; pct(sum(k), length(k)) }
t2row <- function(label, res, term, ihc) tibble(Exposure = label, `IH, % (cases/N)` = ihc,
  Crude = get(res, "Crude", term), `Primary (covariates + site)` = get(res, "Primary", term),
  `Secondary (covariates, no site)` = get(res, "Secondary", term))
table2 <- bind_rows(
  tibble(Exposure = "Single exposures"),
  t2row("Any MVM vs none", r_mvm2, "mvm2n", paste(ih(d1$mvm2n, 1), "vs", ih(d1$mvm2n, 0))),
  t2row("MVM score, per 1 point", r_mvm, "mvm", pct(700, 40700)),
  t2row("Any AI vs none", r_ai, "ain", paste(ih(d1$ain, 1), "vs", ih(d1$ain, 0))),
  t2row("AI score, per 1 compartment", r_aic, "ai_count", pct(700, 40700)),
  tibble(Exposure = "Joint exposure"),
  tibble(Exposure = "Neither (reference)", `IH, % (cases/N)` = ih(d1$joint, "Neither"),
         Crude = "1 REF", `Primary (covariates + site)` = "1 REF", `Secondary (covariates, no site)` = "1 REF"),
  t2row("MVM only", r_joint, "jointMVM only", ih(d1$joint, "MVM only")),
  t2row("AI only", r_joint, "jointAI only", ih(d1$joint, "AI only")),
  t2row("Both MVM and AI", r_joint, "jointBoth MVM and AI", ih(d1$joint, "Both MVM and AI"))) %>%
  mutate(across(everything(), ~ replace_na(.x, "")))
write_csv(table2, file.path(out_dir, paste0("table2_main", sfx, ".csv")))

# Interaction (RERI, AP, S, ratio of RRs) from the joint models
ci <- function(est, g, V, log = FALSE) {
  se <- sqrt(as.numeric(t(g) %*% V %*% g))
  v <- c(est, est - 1.96 * se, est + 1.96 * se); if (log) v <- exp(v)
  sprintf("%.2f (%.2f, %.2f)", v[1], v[2], v[3])
}
jt <- c("jointMVM only", "jointAI only", "jointBoth MVM and AI")
inter <- function(res) {
  b <- res$est$b[match(jt, res$est$term)]; V <- res$V[jt, jt]
  r10 <- exp(b[1]); r01 <- exp(b[2]); r11 <- exp(b[3]); reri <- r11 - r10 - r01 + 1
  c(`RR, MVM only vs neither` = ci(b[1], c(1, 0, 0), V, TRUE),
    `RR, AI only vs neither` = ci(b[2], c(0, 1, 0), V, TRUE),
    `RR, both vs neither` = ci(b[3], c(0, 0, 1), V, TRUE),
    `RERI (additive, 0 = none)` = ci(reri, c(-r10, -r01, r11), V),
    `Attributable proportion, AP (additive, 0 = none)` =
      ci(reri / r11, c(-r10 / r11, -r01 / r11, (r10 + r01 - 1) / r11), V),
    `Synergy index, S (additive, 1 = none)` =
      ci(log((r11 - 1) / (r10 + r01 - 2)), c(-r10, -r01, r11 * (r10 + r01 - 2) / (r11 - 1)) / (r10 + r01 - 2), V, TRUE),
    `Ratio of RRs (multiplicative, 1 = none)` = ci(b[3] - b[1] - b[2], c(-1, -1, 1), V, TRUE))
}
syn <- map_dfc(r_joint, inter)
syn <- tibble(measure = names(inter(r_joint$Crude)), Crude = syn$Crude,
              `Primary: covariates + site + IPW` = syn$Primary,
              `Secondary: covariates + IPW, no site` = syn$Secondary, n_cases = 700, n_obs = 40700)
write_csv(syn, file.path(out_dir, paste0("synergy_reri", sfx, ".csv")))

# Mutual adjustment (same layout as mutual_4exp.csv)
mut <- tibble(
  exposure = c("MVM any vs none", "MVM score, per 1 point", "AI any vs none", "AI score, per 1 compartment"),
  adjusted_for = c("AI any", "AI score", "MVM any", "MVM score"),
  `Crude, mutually adjusted` = c(get(r_mut, "Crude", "mvm2n"), get(r_mutc, "Crude", "mvm"),
                                 get(r_mut, "Crude", "ain"), get(r_mutc, "Crude", "ai_count")),
  `Secondary: + covariates + IPW, no site` = c(get(r_mut, "Secondary", "mvm2n"), get(r_mutc, "Secondary", "mvm"),
                                               get(r_mut, "Secondary", "ain"), get(r_mutc, "Secondary", "ai_count")),
  `Primary: + covariates + site + IPW` = c(get(r_mut, "Primary", "mvm2n"), get(r_mutc, "Primary", "mvm"),
                                           get(r_mut, "Primary", "ain"), get(r_mutc, "Primary", "ai_count")),
  `Primary, single exposure (for comparison)` = c(get(r_mvm2, "Primary", "mvm2n"), get(r_mvm, "Primary", "mvm"),
                                                  get(r_ai, "Primary", "ain"), get(r_aic, "Primary", "ai_count")),
  n_cases = 700, n_obs = 40700)
write_csv(mut, file.path(out_dir, paste0("mutual_4exp", sfx, ".csv")))

# Grade and stage (same layout as grade_stage_separate.csv)
gs_rows <- function(res, x, prefix, label) {
  map_dfr(c("None", "Low", "High"), function(l) {
    k <- d1$ih1[x == l]
    tibble(analysis = label, level = l, ih_pct = pct(sum(k), length(k)), n_cases = sum(k), n_obs = length(k),
           Crude = if (l == "None") "1.00 (ref)" else get(res, "Crude", paste0(prefix, l)),
           `Primary: covariates + site + IPW` = if (l == "None") "1.00 (ref)" else get(res, "Primary", paste0(prefix, l)),
           `Secondary: covariates + IPW, no site` = if (l == "None") "1.00 (ref)" else get(res, "Secondary", paste0(prefix, l)))
  })
}
write_csv(bind_rows(gs_rows(r_g, d1$mvm3, "mvm3", "MVM grade"), gs_rows(r_s, d1$ai3, "ai3", "AI stage")),
          file.path(out_dir, paste0("grade_stage_separate", sfx, ".csv")))

# Maternal or foetal involvement
mf <- map_dfr(levels(d1$aimf), function(l) {
  k <- d1$ih1[d1$aimf == l]
  tibble(`AI location` = l, `IH, % (cases/N)` = pct(sum(k), length(k)),
         Crude = if (l == "No AI") "1 REF" else get(r_mf, "Crude", paste0("aimf", l)),
         Primary = if (l == "No AI") "1 REF" else get(r_mf, "Primary", paste0("aimf", l)),
         Secondary = if (l == "No AI") "1 REF" else get(r_mf, "Secondary", paste0("aimf", l)),
         n_cases = 700, n_obs = 40700)
})
write_csv(mf, file.path(out_dir, paste0("maternal_fetal_ai", sfx, ".csv")))

# Site terms from the primary models (eTable 2)
site_names <- c("5" = "Boston", "10" = "Buffalo", "15" = "New Orleans", "31" = "New York (Columbia)",
                "37" = "Baltimore", "45" = "Virginia", "50" = "Minnesota", "55" = "New York Medical College",
                "60" = "Oregon", "66" = "Pennsylvania", "71" = "Providence", "82" = "Tennessee")
st <- map_dfr(names(site_names), function(s) {
  k <- d1$ih1[d1$site == s]
  rr <- function(res) if (s == "5") "1 REF" else res$Primary$est %>% filter(term == paste0("site", s)) %>% pull(rr_ci)
  tibble(Site = paste(s, site_names[[s]]), `IH cases` = pct(sum(k), length(k)),
         `Site RR (MVM model)` = rr(r_mvm2), `Site RR (AI model)` = rr(r_ai))
})
write_csv(st, file.path(out_dir, paste0("site_table_percentages", sfx, ".csv")))

# ------------------------------------------------------------------------------
# 2. SEQUENTIAL ADJUSTMENT (unweighted cumulative blocks, then sec./primary)
# ------------------------------------------------------------------------------
cat("\n== Sequential\n")
demog <- c("age", "race", "educ", "income", "marital"); health <- c("bmi", "smoking", "dm", "chronic_htn")
repro <- c("parity", "infant_sex")
steps <- list(`1. Crude` = character(0),
  `2. + Demographics (age, race, education, income, marital)` = demog,
  `3. + Health (BMI, smoking, diabetes, chronic HTN)` = c(demog, health),
  `4. + Reproductive (parity, infant sex) = full, no site` = c(demog, health, repro),
  `5. + Study site = full with site` = c(demog, health, repro, "site"))
ex4 <- c(`MVM any vs none` = "mvm2n", `MVM score (continuous)` = "mvm",
         `AI any vs none` = "ain", `AI score (continuous)` = "ai_count")
prim <- list(mvm2n = r_mvm2, mvm = r_mvm, ain = r_ai, ai_count = r_aic)
seq <- imap_dfr(steps, function(cv, s) {
  cat("   ", s, "\n")
  as_tibble(as.list(map_chr(ex4, function(e) fit_pool(e, cv, FALSE, keep = paste0("^", e, "$"))$est$rr_ci)))%>%
    mutate(step = s)
}) %>%
  bind_rows(as_tibble(as.list(map_chr(ex4, ~ get(prim[[.x]], "Secondary", .x)))) %>% mutate(step = "SECONDARY: full, no site, + IPW"),
            as_tibble(as.list(map_chr(ex4, ~ get(prim[[.x]], "Primary", .x)))) %>% mutate(step = "PRIMARY: full + site + IPW")) %>%
  select(step, everything())
write_csv(seq, file.path(out_dir, paste0("seq_cumulative_4exp", sfx, ".csv")))

# ------------------------------------------------------------------------------
# 3. BOSTON VS OTHER 11 CENTRES
# ------------------------------------------------------------------------------
cat("\n== Boston\n")
bos <- lapply(imp, function(d) d[d$site == "5", ])
oth <- lapply(imp, function(d) { d <- d[d$site != "5", ]; d$site <- droplevels(d$site); d })
lab_b <- sprintf("Boston (site 5) [%d cases / %s]", sum(bos[[1]]$ih1), fmt_n(nrow(bos[[1]])))
lab_o <- sprintf("Other 11 sites (+ site) [%d cases / %s]", sum(oth[[1]]$ih1), fmt_n(nrow(oth[[1]])))
bst <- bind_rows(
  as_tibble(as.list(map_chr(ex4, ~ fit_pool(.x, covs, TRUE, data = bos, keep = paste0("^", .x, "$"))$est$rr_ci)))%>% mutate(stratum = lab_b),
  as_tibble(as.list(map_chr(ex4, ~ fit_pool(.x, c(covs, "site"), TRUE, data = oth, keep = paste0("^", .x, "$"))$est$rr_ci)))%>% mutate(stratum = lab_o),
  as_tibble(as.list(map_chr(ex4, ~ get(prim[[.x]], "Primary", .x)))) %>% mutate(stratum = "All sites (primary) [700 cases / 40,700]")) %>%
  select(stratum, everything())
write_csv(bst, file.path(out_dir, paste0("boston_4exp", sfx, ".csv")))
rm(bos, oth); invisible(gc())

# ------------------------------------------------------------------------------
# 4. SUSPECT IH AS IH (ih2), same layout as corrected_outcome_results.csv
# ------------------------------------------------------------------------------
cat("\n== Suspect IH\n")
ex3 <- c(`MVM any vs none` = "mvm2n", `MVM score (continuous)` = "mvm", `AI any vs none` = "ain")
mlab <- c(Crude = "crude", Primary = "adjusted+site+IPW", Secondary = "adjusted+IPW")
ih2 <- imap_dfr(ex3, function(e, lab) imap_dfr(spec, function(s, ml) {
  r <- fit_pool(e, s$covs, s$w, outcome = "ih2", keep = paste0("^", e, "$"))$est
  tibble(analysis = "V2 suspect = present", outcome = "ih2", exposure = lab, model = mlab[[ml]], term = e,
         n_cases = sum(d1$ih2), n_obs = nrow(d1), rr_ci = r$rr_ci, rr = r$rr, lci = r$lci, uci = r$uci, n_imps = r$n_imps)
}))
write_csv(ih2, file.path(out_dir, paste0("ih2", sfx, ".csv")))

# ------------------------------------------------------------------------------
# 5. SITE AS A RANDOM EFFECT (Poisson GLMM with IPW, as 21)
# ------------------------------------------------------------------------------
cat("\n== Site random effect\n")
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
conf_c <- "age_c + bmi_c + race + educ + income + marital + smoking + parity + dm + chronic_htn + infant_sex"
fit_re <- function(e) {
  f <- as.formula(paste("ih1 ~", e, "+", conf_c, "+ (1 | site)"))
  one <- function(d) tryCatch({
    m <- suppressWarnings(glmer(f, data = d, family = poisson(link = "log"), weights = ipw, control = ctrl))
    list(b = fixef(m)[[e]], v = vcov(m)[e, e], s2 = as.numeric(VarCorr(m)$site))
  }, error = function(err) NULL)
  fits <- Filter(Negate(is.null), mclapply(imp, one, mc.cores = n_cores))
  b <- sapply(fits, `[[`, "b"); v <- sapply(fits, `[[`, "v"); m <- length(b)
  se <- sqrt(mean(v) + (1 + 1/m) * var(b))
  list(rr_ci = fmt_ci(exp(mean(b)), exp(mean(b) - 1.96 * se), exp(mean(b) + 1.96 * se)),
       mrr = round(exp(sqrt(2 * mean(sapply(fits, `[[`, "s2"))) * qnorm(0.75)), 2))
}
re_lab <- c(mvm2n = "MVM any vs none", mvm = "MVM score, per 1 point", ain = "AI any vs none", ai_count = "AI score, per 1 compartment")
re <- imap_dfr(re_lab, function(lab, e) {
  cat("   ", lab, "\n"); r <- fit_re(e)
  tibble(exposure = lab, n_cases = 700, n_obs = 40700, Crude = get(prim[[e]], "Crude", e),
         `Primary: site as fixed effect + IPW` = get(prim[[e]], "Primary", e),
         `Site as random effect + IPW` = r$rr_ci,
         `Secondary: no site + IPW` = get(prim[[e]], "Secondary", e), median_rate_ratio_site = r$mrr)
})
write_csv(re, file.path(out_dir, paste0("site_random_effect_models", sfx, ".csv")))

# ------------------------------------------------------------------------------
# 6. E-VALUES (VanderWeele and Ding 2017), primary model
# ------------------------------------------------------------------------------
ev <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
ev_ci <- function(l, u) ifelse(l <= 1 & u >= 1, 1, ifelse(l > 1, ev(l), ev(u)))
pe <- function(res, term, label) {
  x <- res$Primary$est %>% filter(.data$term == !!term)
  tibble(exposure = label, rr_ci = x$rr_ci,
         `E-value, estimate` = round(ev(x$rr), 2), `E-value, CI limit` = round(ev_ci(x$lci, x$uci), 2))
}
evs <- bind_rows(pe(r_joint, "jointBoth MVM and AI", "Both MVM and AI vs neither"),
                 pe(r_joint, "jointMVM only", "MVM only vs neither"),
                 pe(r_joint, "jointAI only", "AI only vs neither"),
                 pe(r_mvm2, "mvm2n", "Any MVM vs none"), pe(r_mvm, "mvm", "MVM score, per point"),
                 pe(r_ai, "ain", "Any AI vs none"), pe(r_aic, "ai_count", "AI score, per compartment"))
write_csv(evs, file.path(out_dir, paste0("evalues", sfx, ".csv")))

write_csv(bind_rows(all_est) %>% filter(!grepl("^site", term)) %>% select(analysis, model, term, rr_ci, rr, lci, uci, n_imps),
          file.path(out_dir, paste0("all_estimates", sfx, ".csv")))

cat("\n==================== TABLE 2 (m = ", M, ") ====================\n", sep = "")
print(table2, width = Inf)
cat("\nInteraction:\n"); print(syn, width = Inf)
cat("\nE-values:\n"); print(evs, width = Inf)
cat("\nSaved _m50 outputs\n")
