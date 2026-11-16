# ==============================================================================
# 70-confounder-drop-check.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Sensitivity check of the confounder set (committee chair's question: which
# covariates are true confounders, and which change the estimates?).
# Exposure: four-level placental injury pattern (neither = reference).
# Uses the existing 50 imputations and stabilised IPCW (49 -> 50), no
# re-imputation. Data prep and fit_pool() as in 51 and 52. Existing results,
# scripts and manuscript files are not changed.
#
# Models (all weighted, pooled with Rubin's rules, 40,700 infants, 700 cases):
#   1  Crude
#   2  Full (12 covariates) -> must reproduce Table 2, else stop
#   3  Drop one at a time: education, income, marital status, smoking,
#      pre-pregnancy BMI, pre-pregnancy diabetes (6 models)
#   4  Drop the socioeconomic block (education, income, marital status)
#   5  Reduced set: drop all 6 (keep age, race, parity, chronic HTN, infant
#      sex, study site)
#   6  As 5, with censoring weights re-estimated from the reduced set
#      (weight model as in 50: denominator = MVM + AI + covariates + site,
#      numerator = intercept only, within each imputation)
#
# Descriptives (observed data, 40,700 examined infants, missing excluded):
#   Cramer's V for education, income and marital status pairs
#   Pre-pregnancy diabetes overall and by category, IH cases among diabetes
#   For the 6 covariates: largest SMD vs neither (as Table 1) and crude RR
#   with IH (unweighted modified Poisson GEE, clustered on mother)
#
# Input:  data/weighted_long_m50y.RDS
# Output: outputs/confounder_drop_check.csv, outputs/confounder_descriptives.csv
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
n_cores <- min(4, detectCores())
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)

# ------------------------------------------------------------------------------
# DATA
# ------------------------------------------------------------------------------
w <- readRDS(file.path(in_dir, "weighted_long_m50y.RDS"))
prep <- function(d) d %>% mutate(
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
  joint   = factor(case_when(mvm2n == 0 & ain == 0 ~ "Neither", mvm2n == 1 & ain == 0 ~ "MVM only",
                             mvm2n == 0 & ain == 1 ~ "AI only", TRUE ~ "Both MVM and AI"),
                   levels = c("Neither", "MVM only", "AI only", "Both MVM and AI"))) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), ~ factor(as.character(.x))))

full_long <- w %>% filter(.imp > 0) %>% prep() %>%
  select(.imp, MOMID, PREGID, CHILDID, ih1, exclusion, ih_obs, ipw, site, age, bmi, race, educ,
         income, marital, smoking, parity, dm, chronic_htn, infant_sex, mvm2, ai, mvm2n, ain, joint)
obs0 <- w %>% filter(.imp == 0, !is.na(ih1), exclusion == 0) %>% prep()   # observed data
rm(w); invisible(gc())

all_imp <- split(full_long, full_long$.imp)                 # 45,268 per imputation (for weights)
imp <- lapply(all_imp, function(d) d %>% filter(!is.na(ih1), exclusion == 0) %>% arrange(MOMID))
d1 <- imp[[1]]; M <- length(imp)
stopifnot(M == 50, nrow(d1) == 40700, sum(d1$ih1) == 700, nrow(all_imp[[1]]) == 45268, nrow(obs0) == 40700)

# ------------------------------------------------------------------------------
# FIT AND POOL (as 51 and 52)
# ------------------------------------------------------------------------------
fit_pool <- function(covs, wcol = "ipw", data = imp) {
  f <- as.formula(paste("ih1 ~", paste(c("joint", covs), collapse = " + ")))
  one <- function(d) {
    d$.w <- d[[wcol]]
    tryCatch({
      m <- geeglm(f, data = d, family = poisson(link = "log"), id = MOMID,
                  corstr = "exchangeable", weights = .w)
      list(b = coef(m), V = vcov(m))
    }, error = function(e) NULL)
  }
  fits <- mclapply(data, one, mc.cores = n_cores)
  for (i in which(sapply(fits, is.null))) fits[[i]] <- one(data[[i]])
  fits <- Filter(Negate(is.null), fits)
  stopifnot(length(fits) == length(data))
  m <- length(fits)
  B <- do.call(rbind, lapply(fits, `[[`, "b"))
  terms <- grep("^joint", colnames(B), value = TRUE)
  B <- B[, terms, drop = FALSE]
  W <- Reduce(`+`, lapply(fits, function(x) x$V[terms, terms, drop = FALSE])) / m
  V <- W + (1 + 1/m) * cov(B)
  b <- colMeans(B); se <- sqrt(diag(V))
  tibble(category = sub("^joint", "", terms), rr = exp(b), lci = exp(b - 1.96 * se), uci = exp(b + 1.96 * se))
}

full_covs <- c("age", "bmi", "race", "educ", "income", "marital", "smoking", "parity",
               "dm", "chronic_htn", "infant_sex", "site")
lab <- c(educ = "education", income = "income", marital = "marital status", smoking = "smoking",
         bmi = "pre-pregnancy BMI", dm = "pre-pregnancy diabetes")
six <- names(lab)
reduced <- setdiff(full_covs, six)

# 2. Full model first: must reproduce Table 2
cat("Fitting full model\n")
r_full <- fit_pool(full_covs)
chk <- sprintf("%.2f (%.2f, %.2f)", r_full$rr, r_full$lci, r_full$uci)
target <- c(`MVM only` = "1.08 (0.89, 1.32)", `AI only` = "1.07 (0.86, 1.32)",
            `Both MVM and AI` = "1.37 (1.05, 1.79)")
print(tibble(category = r_full$category, full = chk, table2 = target[r_full$category]))
if (!all(chk == target[r_full$category])) stop("Full model does not reproduce Table 2. Stopping.")

# Re-estimated weights from the reduced covariate set (model 6), as in 50
rw <- function(d, covs) {
  fd <- as.formula(paste("ih_obs ~ mvm2 + ai +", paste(covs, collapse = " + ")))
  pd <- predict(glm(fd, data = d, family = binomial(link = "logit"), na.action = na.exclude), newdata = d, type = "response")
  pn <- predict(glm(ih_obs ~ 1, data = d, family = binomial(link = "logit")), newdata = d, type = "response")
  as.numeric(pn / pd)
}
# Check the weight code reproduces the stored weights with the full set
stopifnot(isTRUE(all.equal(as.numeric(rw(all_imp[[1]], full_covs)), as.numeric(all_imp[[1]]$ipw), tolerance = 1e-6)))
imp_rw <- lapply(all_imp, function(d) {
  d$ipw_reduced <- rw(d, reduced)
  d %>% filter(!is.na(ih1), exclusion == 0) %>% arrange(MOMID)
})
ws <- imp_rw[[1]]$ipw_reduced
cat(sprintf("Reduced-set weights (imputation 1, examined): mean %.3f, min %.3f, max %.3f\n", mean(ws), min(ws), max(ws)))

models <- c(
  list(`1. Crude` = list(covs = character(0), dropped = "all 12")),
  list(`2. Full (12 covariates)` = list(covs = full_covs, dropped = "none")),
  setNames(lapply(six, function(v) list(covs = setdiff(full_covs, v), dropped = lab[[v]])),
           paste0("3. Drop ", lab[six])),
  list(`4. Drop socioeconomic block` = list(covs = setdiff(full_covs, c("educ", "income", "marital")),
                                            dropped = "education, income, marital status")),
  list(`5. Reduced set (drop all 6)` = list(covs = reduced, dropped = paste(lab[six], collapse = ", "))),
  list(`6. Reduced set, weights re-estimated` = list(covs = reduced, dropped = paste(lab[six], collapse = ", "),
                                                     wcol = "ipw_reduced", data = "rw")))

cache <- file.path(in_dir, "confounder_drop_fits_70.RDS")   # git-ignored, avoids refitting
res_raw <- if (file.exists(cache)) readRDS(cache) else imap_dfr(models, function(s, name) {
  cat("Fitting", name, "\n")
  r <- if (name == "2. Full (12 covariates)") r_full else
    fit_pool(s$covs, wcol = if (is.null(s$wcol)) "ipw" else s$wcol,
             data = if (identical(s$data, "rw")) imp_rw else imp)
  r %>% mutate(model = name, covariates_dropped = s$dropped)
})
saveRDS(res_raw, cache)
res <- res_raw %>%
  left_join(r_full %>% select(category, rr_full = rr), by = "category") %>%
  mutate(pct_change_from_full = 100 * (rr - rr_full) / rr_full,
         flag_10pct = if_else(abs(pct_change_from_full) >= 10 & model != "2. Full (12 covariates)", "YES", ""),
         category = dplyr::recode(category, `MVM only` = "Maternal vascular malperfusion only",
                           `AI only` = "Acute inflammation only", `Both MVM and AI` = "Both"),
         `aRR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", rr, lci, uci),
         n_cases = 700, n_obs = 40700, n_imps = M) %>%
  transmute(model, covariates_dropped, exposure = category, aRR = round(rr, 3), lower_ci = round(lci, 3),
            upper_ci = round(uci, 3), `aRR (95% CI)`, pct_change_from_full = round(pct_change_from_full, 1),
            flag_10pct, n_cases, n_obs, n_imps)
dir.create(out_dir, showWarnings = FALSE)
write_csv(res, file.path(out_dir, "confounder_drop_check.csv"))
print(res %>% select(model, exposure, `aRR (95% CI)`, pct_change_from_full, flag_10pct), n = Inf, width = Inf)

# ------------------------------------------------------------------------------
# DESCRIPTIVES (observed data, missing excluded)
# ------------------------------------------------------------------------------
desc <- list()
add <- function(check, item, value, note = "") desc[[length(desc) + 1]] <<- tibble(check, item, value, note)

cramer <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); tb <- table(droplevels(x[ok]), droplevels(y[ok]))
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
  sqrt(as.numeric(chi) / (sum(tb) * (min(dim(tb)) - 1)))
}
for (p in list(c("educ", "income"), c("educ", "marital"), c("income", "marital")))
  add("Cramer's V", paste(lab[p[1]], "vs", lab[p[2]]), sprintf("%.2f", cramer(obs0[[p[1]]], obs0[[p[2]]])),
      "0 = no association, 1 = complete overlap")

pct <- function(k, n) sprintf("%.1f%% (%s/%s)", 100 * k / n, fmt_n(k), fmt_n(n))
dmv <- obs0$dm == "Yes"
add("Pre-pregnancy diabetes", "Overall", pct(sum(dmv, na.rm = TRUE), sum(!is.na(dmv))))
for (l in levels(obs0$joint)) {
  k <- dmv[obs0$joint == l]; add("Pre-pregnancy diabetes", paste("Category:", l), pct(sum(k, na.rm = TRUE), sum(!is.na(k))))
}
add("Pre-pregnancy diabetes", "IH cases among infants of mothers with diabetes",
    pct(sum(obs0$ih1[dmv %in% TRUE]), sum(dmv, na.rm = TRUE)), "definite IH")
for (l in levels(obs0$joint)) {
  k <- obs0$ih1[dmv %in% TRUE & obs0$joint == l]
  add("Pre-pregnancy diabetes", paste("IH cases with diabetes, category:", l), pct(sum(k), length(k)))
}

smd_cont <- function(x1, x2) abs(mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)) /
  sqrt((var(x1, na.rm = TRUE) + var(x2, na.rm = TRUE)) / 2)
smd_cat <- function(x1, x2) {           # Yang and Dalton (2012), missing excluded
  lv <- levels(x1)
  p1 <- prop.table(table(factor(x1[!is.na(x1)], lv))); p2 <- prop.table(table(factor(x2[!is.na(x2)], lv)))
  if (length(lv) == 2) { a <- p1[2]; b <- p2[2]; return(as.numeric(abs(a - b) / sqrt((a * (1 - a) + b * (1 - b)) / 2))) }
  dd <- (p1 - p2)[-1]; q1 <- p1[-1]; q2 <- p2[-1]
  S <- (-(outer(q1, q1) + outer(q2, q2)) / 2); diag(S) <- (q1 * (1 - q1) + q2 * (1 - q2)) / 2
  as.numeric(sqrt(t(dd) %*% solve(S) %*% dd))
}
crude_rr <- function(var) {
  d <- obs0 %>% filter(!is.na(.data[[var]])) %>% arrange(MOMID)
  x <- if (var == "bmi") "I(bmi / 5)" else var
  m <- geeglm(as.formula(paste("ih1 ~", x)), data = d, family = poisson(link = "log"), id = MOMID, corstr = "exchangeable")
  s <- summary(m)$coefficients[-1, , drop = FALSE]
  tibble(term = rownames(s), rr_ci = sprintf("%.2f (%.2f, %.2f)", exp(s[, 1]), exp(s[, 1] - 1.96 * s[, 2]), exp(s[, 1] + 1.96 * s[, 2])))
}
for (v in six) {
  f <- if (v == "bmi") smd_cont else smd_cat
  s <- map_dbl(c("MVM only", "AI only", "Both MVM and AI"), ~ f(obs0[[v]][obs0$joint == .x], obs0[[v]][obs0$joint == "Neither"]))
  add("Largest SMD vs neither", lab[[v]], sprintf("%.2f", max(s)),
      sprintf("MVM only %.2f, AI only %.2f, both %.2f", s[1], s[2], s[3]))
  ref <- if (v == "bmi") "per 5 kg/m2" else paste("reference:", levels(obs0[[v]])[1])
  cr <- crude_rr(v)
  for (i in seq_len(nrow(cr))) add("Crude RR with IH (unweighted)", paste0(lab[[v]], ": ", cr$term[i]), cr$rr_ci[i], ref)
}
desc <- bind_rows(desc)
write_csv(desc, file.path(out_dir, "confounder_descriptives.csv"))
print(desc, n = Inf, width = Inf)
cat("\nDone\n")
