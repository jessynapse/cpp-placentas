# ==============================================================================
# 71-reduced-covariates.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Reduced confounder set (committee chair's question about true confounders).
# Reduced set (9): maternal age, race and ethnicity, education, income,
# marital status, smoking, parity, infant sex, study site. Removed from the
# primary 12: pre-pregnancy BMI, pre-pregnancy diabetes, chronic hypertension.
#
# Everything else as the current pipeline (49 -> 50 -> 51, 52, 54, 55):
# 40,700 examined infants, 700 definite IH, same 50 imputations (no
# re-imputation), modified Poisson GEE, exchangeable, clustered on MOMID,
# Rubin's rules with the pooled covariance, delta-method CIs.
# Stabilised IPCW re-estimated within each imputation on all 45,268 pregnancies
# (fetal deaths, infant deaths and survival without examination = censored),
# denominator = mvm2 + ai + covariates, numerator = intercept only (as 50).
# Crude models are unweighted with no covariates (as 51), so they are the same
# for both sets. The 12-covariate results are refitted here with the stored
# weights and checked against v4.6. If Table 2 does not reproduce, it stops.
#
# Tables rerun (v4.6 numbering):
#   Table 2   four-category pattern, crude and adjusted, RERI, AP, ratio of RRs
#   eTable 2  each pattern within strata of the other (from the joint model)
#   eTable 3  each pattern any vs none, MVM grade, AI stage
#   eTable 4  Boston site (no site term) vs other 11 sites (site term)
#   eTable 5  definite IH vs definite or suspect IH (ih2)
#   eTable 6  E-values and the maternal angioma bound (RR_UD = 1.64)
# Intermediate models (four-category pattern only, weights re-estimated with
# the same covariates): 12 minus chronic hypertension, and 12 minus BMI and
# diabetes.
#
# Input:  data/weighted_long_m50y.RDS
# Output: outputs/reduced9_comparison.csv, outputs/table2_reduced9.csv,
#         outputs/etable2_reduced9.csv ... outputs/etable6_reduced9.csv
# Cache:  data/reduced9_fits_71.RDS (git-ignored, pooled fits)
# Test:   R71_TEST=3 Rscript ... uses 3 imputations and writes nothing to outputs
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
test_n  <- as.integer(Sys.getenv("R71_TEST", "0"))
if (test_n > 0) out_dir <- tempdir()
n_cores <- min(4, detectCores())
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)
pct   <- function(k, n) sprintf("%.1f%% (%s/%s)", 100 * k / n, fmt_n(k), fmt_n(n))
f2    <- function(r, l, u) sprintf("%.2f (%.2f, %.2f)", r, l, u)
REF   <- "1.00 (Reference)"
fx    <- function(x) f2(x[[1]], x[[2]], x[[3]])

# ------------------------------------------------------------------------------
# DATA (factor coding as 51 and 52)
# ------------------------------------------------------------------------------
w <- readRDS(file.path(in_dir, "weighted_long_m50y.RDS"))
if (test_n > 0) w <- w[w$.imp <= test_n, ]
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
  mvm3    = factor(as.character(mvm3), levels = c("None", "Low", "High")),
  ai3     = factor(as.character(ai3),  levels = c("None", "Low", "High")),
  joint   = factor(case_when(mvm2n == 0 & ain == 0 ~ "Neither", mvm2n == 1 & ain == 0 ~ "MVM only",
                             mvm2n == 0 & ain == 1 ~ "AI only", TRUE ~ "Both MVM and AI"),
                   levels = c("Neither", "MVM only", "AI only", "Both MVM and AI"))) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), ~ factor(as.character(.x))))

all_imp <- w %>% filter(.imp > 0) %>% prep() %>%
  select(.imp, MOMID, ih1, ih2, exclusion, ih_obs, ipw, site, age, bmi, race, educ, income, marital,
         smoking, parity, dm, chronic_htn, infant_sex, mvm2, ai, mvm2n, ain, mvm3, ai3, joint) %>%
  split(.$.imp)
rm(w); invisible(gc())
M <- length(all_imp)
stopifnot(test_n > 0 || M == 50, all(sapply(all_imp, nrow) == 45268))

full12 <- c("age", "bmi", "race", "educ", "income", "marital", "smoking", "parity",
            "dm", "chronic_htn", "infant_sex", "site")
red9   <- setdiff(full12, c("bmi", "dm", "chronic_htn"))
noHTN  <- setdiff(full12, "chronic_htn")
noBMDM <- setdiff(full12, c("bmi", "dm"))
stopifnot(length(red9) == 9)

# Stabilised IPCW, as 50
rw <- function(d, covs) {
  fd <- as.formula(paste("ih_obs ~ mvm2 + ai +", paste(covs, collapse = " + ")))
  pd <- predict(glm(fd, data = d, family = binomial(link = "logit"), na.action = na.exclude), newdata = d, type = "response")
  pn <- predict(glm(ih_obs ~ 1, data = d, family = binomial(link = "logit")), newdata = d, type = "response")
  as.numeric(pn / pd)
}
stopifnot(isTRUE(all.equal(rw(all_imp[[1]], full12), as.numeric(all_imp[[1]]$ipw), tolerance = 1e-6)))
cat("Stored weights reproduce from the 12-covariate weight model\n")

imp <- lapply(all_imp, function(d) {
  d$ipw9 <- rw(d, red9); d$ipw_noHTN <- rw(d, noHTN); d$ipw_noBMDM <- rw(d, noBMDM)
  d %>% filter(!is.na(ih1), exclusion == 0) %>% arrange(MOMID)
})
rm(all_imp); invisible(gc())
d1 <- imp[[1]]
stopifnot(nrow(d1) == 40700, sum(d1$ih1) == 700, sum(d1$ih2) == 930)
wsum <- map_dfr(c("ipw", "ipw9", "ipw_noHTN", "ipw_noBMDM"), ~ tibble(weights = .x,
  mean = mean(d1[[.x]]), min = min(d1[[.x]]), max = max(d1[[.x]])))
print(wsum)
bos <- lapply(imp, function(d) d[d$site == "5", ])
oth <- lapply(imp, function(d) { d <- d[d$site != "5", ]; d$site <- droplevels(d$site); d })
stopifnot(sum(bos[[1]]$ih1) == 290, sum(oth[[1]]$ih1) == 410)

# ------------------------------------------------------------------------------
# FIT AND POOL (as 51), with a cache of pooled fits
# ------------------------------------------------------------------------------
cache_f <- file.path(in_dir, "reduced9_fits_71.RDS")
cache <- if (test_n == 0 && file.exists(cache_f)) readRDS(cache_f) else list()

fit_pool <- function(key, expo, covs, wcol = NULL, outcome = "ih1", data = imp) {
  if (!is.null(cache[[key]])) return(cache[[key]])
  cat("Fitting", key, "\n")
  f <- as.formula(paste(outcome, "~", paste(c(expo, covs), collapse = " + ")))
  one <- function(d) {
    d$.w <- if (is.null(wcol)) rep(1, nrow(d)) else d[[wcol]]
    tryCatch({
      m <- geeglm(f, data = d, family = poisson(link = "log"), id = MOMID, corstr = "exchangeable", weights = .w)
      list(b = coef(m), V = vcov(m))
    }, error = function(e) NULL)
  }
  fits <- mclapply(data, one, mc.cores = n_cores)
  for (i in which(sapply(fits, is.null))) fits[[i]] <- one(data[[i]])
  fits <- Filter(Negate(is.null), fits)
  stopifnot(length(fits) == length(data))
  m <- length(fits)
  B <- do.call(rbind, lapply(fits, `[[`, "b"))
  terms <- grep(paste0("^", expo), colnames(B), value = TRUE)
  B <- B[, terms, drop = FALSE]
  W <- Reduce(`+`, lapply(fits, function(x) x$V[terms, terms, drop = FALSE])) / m
  V <- W + (1 + 1/m) * cov(B)
  res <- list(b = colMeans(B), V = V, n_cases = sum(data[[1]][[outcome]]), n_obs = nrow(data[[1]]), n_imps = m)
  cache[[key]] <<- res
  if (test_n == 0) saveRDS(cache, cache_f)
  res
}
# Ratio estimates from a pooled fit: g = contrast on the log scale (exp(g'b))
lin <- function(fit, g) { e <- sum(g * fit$b); s <- sqrt(as.numeric(t(g) %*% fit$V %*% g))
  c(rr = exp(e), lci = exp(e - 1.96 * s), uci = exp(e + 1.96 * s)) }
term_rr <- function(fit, term) lin(fit, as.numeric(names(fit$b) == term))
inter <- function(fit) {           # RERI and AP on the additive scale, ratio of RRs on the log scale
  r <- exp(fit$b[c("jointMVM only", "jointAI only", "jointBoth MVM and AI")])
  V <- fit$V[names(r), names(r)]
  dm <- function(est, g) { s <- sqrt(as.numeric(t(g) %*% V %*% g)); c(rr = est, lci = est - 1.96 * s, uci = est + 1.96 * s) }
  reri <- r[3] - r[1] - r[2] + 1
  list(RERI = dm(reri, c(-r[1], -r[2], r[3])),
       AP   = dm(reri / r[3], c(-r[1] / r[3], -r[2] / r[3], (r[1] + r[2] - 1) / r[3])),
       `Ratio of RRs` = { g <- c(-1, -1, 1); e <- sum(g * log(r)); s <- sqrt(as.numeric(t(g) %*% V %*% g))
                          c(rr = exp(e), lci = exp(e - 1.96 * s), uci = exp(e + 1.96 * s)) })
}
jt <- c(`Maternal vascular malperfusion only` = "jointMVM only", `Acute inflammation only` = "jointAI only",
        Both = "jointBoth MVM and AI")
gw <- function(fit, a, b) { g <- setNames(rep(0, length(fit$b)), names(fit$b)); g[a] <- 1; if (!is.null(b)) g[b] <- -1; lin(fit, g) }
strata <- list(
  `MVM, among placentas without AI` = c("jointMVM only"),
  `MVM, among placentas with AI`    = c("jointBoth MVM and AI", "jointAI only"),
  `AI, among placentas without MVM` = c("jointAI only"),
  `AI, among placentas with MVM`    = c("jointBoth MVM and AI", "jointMVM only"))
ev    <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
ev_ci <- function(l, u) ifelse(l <= 1 & u >= 1, 1, ifelse(l > 1, ev(l), ev(u)))
bound <- function(rr, ud = 1.64) rr * (ud - 1) / (ud - rr)

# ------------------------------------------------------------------------------
# STEP 0: the 12-covariate pipeline must reproduce Table 2 (v4.6)
# ------------------------------------------------------------------------------
crude_j <- fit_pool("crude_joint", "joint", character(0))
j12 <- fit_pool("12_joint", "joint", full12, "ipw")
got <- c(map_chr(jt, ~ fx(term_rr(j12, .x))),
         map_chr(jt, ~ fx(term_rr(crude_j, .x))),
         fx(inter(j12)$RERI))
want <- c("1.08 (0.89, 1.32)", "1.07 (0.86, 1.32)", "1.37 (1.05, 1.79)",
          "1.18 (0.98, 1.42)", "1.19 (0.96, 1.47)", "1.55 (1.20, 1.99)", "0.22 (-0.20, 0.64)")
print(tibble(check = c(paste("adjusted", names(jt)), paste("crude", names(jt)), "adjusted RERI"), got, want))
if (!all(got == want)) {
  if (test_n > 0) message("TEST MODE: Table 2 not reproduced with ", M, " imputations (expected)") else
    stop("The 12-covariate pipeline does not reproduce Table 2. Stopping before any reduced-set model.")
} else cat("Table 2 reproduced. Continuing.\n")

# ------------------------------------------------------------------------------
# FITS: 12 (stored weights) and 9 (re-estimated weights)
# ------------------------------------------------------------------------------
sets <- list(`12` = list(covs = full12, w = "ipw"), `9` = list(covs = red9, w = "ipw9"))
crude <- list(joint = crude_j,
              mvm2n = fit_pool("crude_mvm2n", "mvm2n", character(0)),
              ain   = fit_pool("crude_ain",   "ain",   character(0)),
              mvm3  = fit_pool("crude_mvm3",  "mvm3",  character(0)),
              ai3   = fit_pool("crude_ai3",   "ai3",   character(0)))
fits <- imap(sets, function(s, k) list(
  joint = if (k == "12") j12 else fit_pool(paste0(k, "_joint"), "joint", s$covs, s$w),
  mvm2n = fit_pool(paste0(k, "_mvm2n"), "mvm2n", s$covs, s$w),
  ain   = fit_pool(paste0(k, "_ain"),   "ain",   s$covs, s$w),
  mvm3  = fit_pool(paste0(k, "_mvm3"),  "mvm3",  s$covs, s$w),
  ai3   = fit_pool(paste0(k, "_ai3"),   "ai3",   s$covs, s$w),
  bos   = fit_pool(paste0(k, "_bos"),   "joint", setdiff(s$covs, "site"), s$w, data = bos),
  oth   = fit_pool(paste0(k, "_oth"),   "joint", s$covs, s$w, data = oth),
  ih2   = fit_pool(paste0(k, "_ih2"),   "joint", s$covs, s$w, outcome = "ih2")))
mid <- list(`12 minus chronic hypertension` = fit_pool("noHTN_joint", "joint", noHTN, "ipw_noHTN"),
            `12 minus pre-pregnancy BMI and diabetes` = fit_pool("noBMDM_joint", "joint", noBMDM, "ipw_noBMDM"))

# ------------------------------------------------------------------------------
# LONG RESULTS: one row per table, exposure and covariate set
# ------------------------------------------------------------------------------
L <- list()
put <- function(table, exposure, set, x, n_cases, n_obs, measure = "RR")
  L[[length(L) + 1]] <<- tibble(table, exposure, set, measure, est = x[[1]], lci = x[[2]], uci = x[[3]], n_cases, n_obs)
for (k in names(fits)) {
  F <- fits[[k]]
  for (e in names(jt)) put("Table 2", e, k, term_rr(F$joint, jt[[e]]), 700, 40700)
  for (e in names(inter(F$joint))) put("Table 2", e, k, inter(F$joint)[[e]], 700, 40700,
                                       if (e == "Ratio of RRs") "RR" else "additive")
  for (e in names(strata)) put("eTable 2", e, k, gw(F$joint, strata[[e]][1], if (length(strata[[e]]) > 1) strata[[e]][2]), 700, 40700)
  put("eTable 3", "Any MVM vs none", k, term_rr(F$mvm2n, "mvm2n"), 700, 40700)
  put("eTable 3", "MVM grade low vs none", k, term_rr(F$mvm3, "mvm3Low"), 700, 40700)
  put("eTable 3", "MVM grade high vs none", k, term_rr(F$mvm3, "mvm3High"), 700, 40700)
  put("eTable 3", "Any AI vs none", k, term_rr(F$ain, "ain"), 700, 40700)
  put("eTable 3", "AI stage low vs none", k, term_rr(F$ai3, "ai3Low"), 700, 40700)
  put("eTable 3", "AI stage high vs none", k, term_rr(F$ai3, "ai3High"), 700, 40700)
  for (e in names(jt)) put("eTable 4", paste("Boston:", e), k, term_rr(F$bos, jt[[e]]), F$bos$n_cases, F$bos$n_obs)
  for (e in names(jt)) put("eTable 4", paste("Other 11 sites:", e), k, term_rr(F$oth, jt[[e]]), F$oth$n_cases, F$oth$n_obs)
  for (e in names(jt)) put("eTable 5", paste("Definite or suspect IH:", e), k, term_rr(F$ih2, jt[[e]]), 930, 40700)
  for (e in names(jt)) { x <- term_rr(F$joint, jt[[e]])
    put("eTable 6", paste("E-value:", e), k, c(ev(x[["rr"]]), ev_ci(x[["lci"]], x[["uci"]]), NA), 700, 40700, "E-value") }
  xb <- term_rr(F$joint, jt[["Both"]])
  put("eTable 6", "Maternal angioma bound (RR_UD 1.64)", k, c(bound(xb[["rr"]]), bound(xb[["lci"]]), NA), 700, 40700, "bound")
}
for (k in names(mid)) for (e in names(jt)) put("Table 2 (intermediate)", e, k, term_rr(mid[[k]], jt[[e]]), 700, 40700)
long <- bind_rows(L)

# ------------------------------------------------------------------------------
# COMPARISON: 12 vs 9, and 12 vs each intermediate
# ------------------------------------------------------------------------------
fmt_row <- function(m, e, l, u) case_when(m == "E-value" ~ sprintf("%.2f (CI bound %.2f)", e, l),
                                          m == "bound" ~ sprintf("%.2f (lower CI %.2f)", e, l),
                                          TRUE ~ f2(e, l, u))
base <- long %>% filter(set == "12") %>% select(table, exposure, measure, est12 = est, lci12 = lci, uci12 = uci)
cmp_main <- long %>% filter(set == "9") %>% mutate(model = "9 covariates (drop BMI, diabetes, chronic hypertension)") %>%
  left_join(base, by = c("table", "exposure", "measure"))
cmp_mid <- long %>% filter(table == "Table 2 (intermediate)") %>% mutate(model = set) %>%
  left_join(base %>% filter(table == "Table 2") %>% select(-table), by = c("exposure", "measure"))
comparison <- bind_rows(cmp_main, cmp_mid) %>%
  mutate(pct_change = ifelse(measure == "additive", NA, 100 * (est - est12) / est12),
         flag_10pct = if_else(!is.na(pct_change) & abs(pct_change) >= 10, "YES", ""),
         `12 covariates` = fmt_row(measure, est12, lci12, uci12),
         `Comparison model` = fmt_row(measure, est, lci, uci)) %>%
  transmute(model, table, exposure, measure, `12 covariates`, `Comparison model`,
            aRR_12 = round(est12, 3), lci_12 = round(lci12, 3), uci_12 = round(uci12, 3),
            aRR_comp = round(est, 3), lci_comp = round(lci, 3), uci_comp = round(uci, 3),
            pct_change = round(pct_change, 1), flag_10pct, n_cases, n_obs, n_imps = M)
write_csv(comparison, file.path(out_dir, "reduced9_comparison.csv"))

# Check the 12-covariate column against v4.6 (every adjusted cell)
v46 <- c("1.08 (0.89, 1.32)", "1.07 (0.86, 1.32)", "1.37 (1.05, 1.79)", "0.22 (-0.20, 0.64)", "0.16 (-0.12, 0.44)",
         "1.19 (0.84, 1.68)", "1.08 (0.89, 1.32)", "1.29 (0.95, 1.75)", "1.07 (0.86, 1.32)", "1.27 (0.96, 1.67)",
         "1.14 (0.96, 1.35)", "1.15 (0.96, 1.37)", "1.09 (0.76, 1.54)", "1.14 (0.96, 1.35)", "1.26 (1.02, 1.55)",
         "1.00 (0.79, 1.28)", "0.91 (0.67, 1.24)", "1.05 (0.77, 1.43)", "1.65 (1.16, 2.34)", "1.17 (0.90, 1.52)",
         "1.08 (0.81, 1.44)", "1.12 (0.75, 1.66)", "1.09 (0.92, 1.29)", "1.05 (0.88, 1.26)", "1.29 (1.03, 1.62)",
         "1.39 (CI bound 1.00)", "1.33 (CI bound 1.00)", "2.09 (CI bound 1.29)", "3.28 (lower CI 1.15)")
chk12 <- tibble(row = paste(cmp_main$table, cmp_main$exposure), current_v46 = v46,
                refit_12 = comparison$`12 covariates`[seq_len(nrow(cmp_main))]) %>%
  mutate(match = current_v46 == refit_12)
print(chk12, n = Inf, width = Inf)
if (test_n == 0 && !all(chk12$match)) warning("Some 12-covariate eTable cells do not match v4.6. See the check above.")

# ------------------------------------------------------------------------------
# TABLES WITH THE 9-COVARIATE SET (layout as v4.6)
# ------------------------------------------------------------------------------
g9 <- function(tab, e) { x <- long %>% filter(set == "9", table == tab, exposure == e); f2(x$est, x$lci, x$uci) }
cr <- function(fit, term) fx(term_rr(fit, term))
cnt <- function(x, y = d1$ih1) pct(sum(y[x]), sum(x))
jl <- c(Neither = "Neither", `Maternal vascular malperfusion only` = "MVM only",
        `Acute inflammation only` = "AI only", Both = "Both MVM and AI")

t2 <- bind_rows(
  tibble(Exposure = names(jl), `IH, % (cases/N)` = map_chr(jl, ~ cnt(d1$joint == .x)),
         `Unadjusted RR (95% CI)` = c(REF, map_chr(jt, ~ cr(crude_j, .x))),
         `aRR (95% CI), 9 covariates` = c(REF, map_chr(names(jt), ~ g9("Table 2", .x)))),
  tibble(Exposure = c("RERI (additive, 0 = none)", "AP (additive, 0 = none)", "Ratio of RRs (multiplicative, 1 = none)"),
         `IH, % (cases/N)` = "",
         `Unadjusted RR (95% CI)` = map_chr(inter(crude_j), ~ fx(.x)),
         `aRR (95% CI), 9 covariates` = map_chr(c("RERI", "AP", "Ratio of RRs"), ~ g9("Table 2", .x))))

st_n <- list(`MVM, among placentas without AI` = c("Neither", "MVM only"), `MVM, among placentas with AI` = c("AI only", "Both MVM and AI"),
             `AI, among placentas without MVM` = c("Neither", "AI only"), `AI, among placentas with MVM` = c("MVM only", "Both MVM and AI"))
et2 <- imap_dfr(strata, function(s, lab) tibble(
  Stratum = lab, Level = c("None", "Any"), `IH, % (cases/N)` = map_chr(st_n[[lab]], ~ cnt(d1$joint == .x)),
  `Unadjusted RR (95% CI)` = c(REF, fx(gw(crude_j, s[1], if (length(s) > 1) s[2])))
  , `aRR (95% CI), 9 covariates` = c(REF, g9("eTable 2", lab))))

et3 <- bind_rows(
  tibble(Exposure = "MVM, any vs none", Level = c("None", "Any"), `IH, % (cases/N)` = c(cnt(d1$mvm2n == 0), cnt(d1$mvm2n == 1)),
         `Unadjusted RR (95% CI)` = c(REF, cr(crude$mvm2n, "mvm2n")), `aRR (95% CI), 9 covariates` = c(REF, g9("eTable 3", "Any MVM vs none"))),
  tibble(Exposure = "MVM grade", Level = c("None", "Low (score 1-2)", "High (score 3 or more)"),
         `IH, % (cases/N)` = map_chr(c("None", "Low", "High"), ~ cnt(d1$mvm3 == .x)),
         `Unadjusted RR (95% CI)` = c(REF, cr(crude$mvm3, "mvm3Low"), cr(crude$mvm3, "mvm3High")),
         `aRR (95% CI), 9 covariates` = c(REF, g9("eTable 3", "MVM grade low vs none"), g9("eTable 3", "MVM grade high vs none"))),
  tibble(Exposure = "AI, any vs none", Level = c("None", "Any"), `IH, % (cases/N)` = c(cnt(d1$ain == 0), cnt(d1$ain == 1)),
         `Unadjusted RR (95% CI)` = c(REF, cr(crude$ain, "ain")), `aRR (95% CI), 9 covariates` = c(REF, g9("eTable 3", "Any AI vs none"))),
  tibble(Exposure = "AI stage", Level = c("None", "Low", "High"),
         `IH, % (cases/N)` = map_chr(c("None", "Low", "High"), ~ cnt(d1$ai3 == .x)),
         `Unadjusted RR (95% CI)` = c(REF, cr(crude$ai3, "ai3Low"), cr(crude$ai3, "ai3High")),
         `aRR (95% CI), 9 covariates` = c(REF, g9("eTable 3", "AI stage low vs none"), g9("eTable 3", "AI stage high vs none"))))

bd <- bos[[1]]; od <- oth[[1]]
et4 <- tibble(Exposure = c("All infants", names(jl)),
  `Boston: IH, % (cases/N)` = c(pct(sum(bd$ih1), nrow(bd)), map_chr(jl, ~ pct(sum(bd$ih1[bd$joint == .x]), sum(bd$joint == .x)))),
  `Boston: aRR (95% CI), 9 covariates` = c("", REF, map_chr(names(jt), ~ g9("eTable 4", paste("Boston:", .x)))),
  `Other 11 sites: IH, % (cases/N)` = c(pct(sum(od$ih1), nrow(od)), map_chr(jl, ~ pct(sum(od$ih1[od$joint == .x]), sum(od$joint == .x)))),
  `Other 11 sites: aRR (95% CI), 9 covariates` = c("", REF, map_chr(names(jt), ~ g9("eTable 4", paste("Other 11 sites:", .x)))))

et5 <- tibble(Exposure = c("All infants", names(jl)),
  `Definite IH: IH, % (cases/N)` = c(pct(700, 40700), map_chr(jl, ~ cnt(d1$joint == .x))),
  `Definite IH: aRR (95% CI), 9 covariates` = c("", REF, map_chr(names(jt), ~ g9("Table 2", .x))),
  `Definite or suspect IH: IH, % (cases/N)` = c(pct(930, 40700), map_chr(jl, ~ cnt(d1$joint == .x, d1$ih2))),
  `Definite or suspect IH: aRR (95% CI), 9 covariates` = c("", REF, map_chr(names(jt), ~ g9("eTable 5", paste("Definite or suspect IH:", .x)))))

e9 <- long %>% filter(set == "9", table == "eTable 6")
et6 <- tibble(Exposure = paste(names(jt), "vs neither"),
  `aRR (95% CI), 9 covariates` = map_chr(names(jt), ~ g9("Table 2", .x)),
  `E-value, point estimate` = sprintf("%.2f", e9$est[1:3]), `E-value, CI bound` = sprintf("%.2f", e9$lci[1:3])) %>%
  bind_rows(tibble(Exposure = "Maternal history of angioma (RR_UD 1.64): RR_EU to explain away the aRR for both",
                   `aRR (95% CI), 9 covariates` = "", `E-value, point estimate` = sprintf("%.2f", e9$est[4]),
                   `E-value, CI bound` = sprintf("%.2f", e9$lci[4])))

note <- "Adjusted for maternal age, race and ethnicity, education, annual family income, marital status, smoking, parity, infant sex and study site, with the same 50 imputations and stabilised censoring weights re-estimated from both patterns and these 9 covariates. 700 definite IH cases among 40,700 infants unless stated."
wr <- function(x, f) write_csv(x %>% mutate(note = c(note, rep("", nrow(x) - 1))), file.path(out_dir, f))
wr(t2, "table2_reduced9.csv"); wr(et2, "etable2_reduced9.csv"); wr(et3, "etable3_reduced9.csv")
wr(et4, "etable4_reduced9.csv"); wr(et5, "etable5_reduced9.csv"); wr(et6, "etable6_reduced9.csv")
write_csv(wsum, file.path(out_dir, "reduced9_weights.csv"))

print(comparison %>% select(model, table, exposure, `12 covariates`, `Comparison model`, pct_change, flag_10pct), n = Inf, width = Inf)
for (x in list(t2, et2, et3, et4, et5, et6)) print(x, n = Inf, width = Inf)
cat("\nDone\n")
