# ==============================================================================
# 63-boston-heterogeneity-clean.R
# COPY of 53-boston-heterogeneity.R for the rerun after 58 (implausible values set to missing).
# Only change: file names (_m50y -> _m50c, analytic_sample_corrected ->
# analytic_sample_clean, and c_ prefixes on CSV outputs), so the earlier
# results are kept for comparison. Code otherwise identical.
# ==============================================================================

# ==============================================================================
# 53-boston-heterogeneity.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Does the joint MVM x AI association differ between Boston and the other
# 11 sites? Follows 52 (Rerun 2), which found interaction in Boston only.
# Same data prep, weights and fit_pool() as 51 and 52. Primary analysis
# unchanged. No p-values: heterogeneity is reported as ratios or
# differences with 95% CIs.
#
# A. Stratified (as 52): the two strata share no mothers, so their
#    estimates are independent and variances add.
#      - Ratio of stratum aRRs (Boston / other sites) for each joint category
#      - Ratio of the two ratios of RRs (multiplicative interaction)
#      - Difference in RERI (additive interaction)
#    Boston has no site term. Other sites: primary with the 11-level site
#    term, secondary without. Boston primary = Boston secondary.
#
# B. Single model, all 40,700: joint + product terms joint x Boston
#    (jb_*) + covariates. Covariate effects are shared across strata.
#      - Primary: + site (Boston main effect is the site 5 term)
#      - Secondary: + Boston indicator, no site
#      - Crude: + Boston indicator only
#    exp(jb_k) = ratio of the category aRR in Boston vs other sites.
#    exp(jb_both - jb_mvm - jb_ai) = ratio of the two ratios of RRs.
#
# Input:  data/weighted_long_m50c.RDS
# Output: outputs/boston_heterogeneity_stratified.csv,
#         outputs/boston_heterogeneity_model.csv
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
n_cores <- min(4, detectCores())   # set to 1 on Windows

fmt_ci <- function(rr, l, u) sprintf("%.2f (%.2f, %.2f)", rr, l, u)

# ------------------------------------------------------------------------------
# DATA (as 51 and 52)
# ------------------------------------------------------------------------------
w <- readRDS(file.path(in_dir, "weighted_long_m50c.RDS"))
imp <- w %>%
  filter(.imp > 0, !is.na(ih1), exclusion == 0) %>%
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
    joint   = factor(case_when(mvm2n == 0 & ain == 0 ~ "Neither", mvm2n == 1 & ain == 0 ~ "MVM only",
                               mvm2n == 0 & ain == 1 ~ "AI only", TRUE ~ "Both MVM and AI"),
                     levels = c("Neither", "MVM only", "AI only", "Both MVM and AI")),
    boston  = as.numeric(site == "5"),
    jb_mvm  = as.numeric(joint == "MVM only") * boston,
    jb_ai   = as.numeric(joint == "AI only") * boston,
    jb_both = as.numeric(joint == "Both MVM and AI") * boston
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), ~ factor(as.character(.x)))) %>%
  select(.imp, MOMID, ih1, ipw, site, boston, age, bmi, race, educ, income, marital,
         smoking, parity, dm, chronic_htn, infant_sex, joint, jb_mvm, jb_ai, jb_both) %>%
  arrange(MOMID) %>%
  split(.$.imp)
rm(w); invisible(gc())
d1 <- imp[[1]]
M <- length(imp)
stopifnot(nrow(d1) == 40700, sum(d1$ih1) == 700, M == 50, sum(d1$ih1[d1$boston == 1]) == 290)

# ------------------------------------------------------------------------------
# FIT AND POOL (copied from 51)
# ------------------------------------------------------------------------------
covs <- c("age", "bmi", "race", "educ", "income", "marital", "smoking", "parity",
          "dm", "chronic_htn", "infant_sex")

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

# Estimate and variance of a function of the joint coefficients (delta method)
jt <- c("jointMVM only", "jointAI only", "jointBoth MVM and AI")
measures <- function(res) {
  b <- res$est$b[match(jt, res$est$term)]; V <- res$V[jt, jt]
  r10 <- exp(b[1]); r01 <- exp(b[2]); r11 <- exp(b[3])
  q <- function(est, g) list(est = est, var = as.numeric(t(g) %*% V %*% g))
  list(`MVM only` = q(b[1], c(1, 0, 0)),                     # log scale
       `AI only` = q(b[2], c(0, 1, 0)),                      # log scale
       `Both MVM and AI` = q(b[3], c(0, 0, 1)),              # log scale
       ratio = q(b[3] - b[1] - b[2], c(-1, -1, 1)),          # log scale
       reri = q(r11 - r10 - r01 + 1, c(-r10, -r01, r11)))    # natural scale
}
diff_ci <- function(a, b, log = TRUE) {
  est <- a$est - b$est; se <- sqrt(a$var + b$var)
  v <- c(est, est - 1.96 * se, est + 1.96 * se); if (log) v <- exp(v)
  sprintf("%.2f (%.2f, %.2f)", v[1], v[2], v[3])
}
one_ci <- function(a, log = TRUE) {
  se <- sqrt(a$var); v <- c(a$est, a$est - 1.96 * se, a$est + 1.96 * se); if (log) v <- exp(v)
  sprintf("%.2f (%.2f, %.2f)", v[1], v[2], v[3])
}

dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# A. STRATIFIED COMPARISON
# ------------------------------------------------------------------------------
cat("\n== A. Stratified\n")
bos <- lapply(imp, function(d) d[d$site == "5", ])
oth <- lapply(imp, function(d) { d <- d[d$site != "5", ]; d$site <- droplevels(d$site); d })
fitj <- function(data, cv, wt) fit_pool("joint", cv, wt, data = data, keep = "^joint")

cat("   Boston\n")
mb <- list(Crude = measures(fitj(bos, character(0), FALSE)),
           Primary = measures(fitj(bos, covs, TRUE)))
mb$Secondary <- mb$Primary   # Boston has no site term in either model
cat("   Other 11 sites\n")
mo <- list(Crude = measures(fitj(oth, character(0), FALSE)),
           Primary = measures(fitj(oth, c(covs, "site"), TRUE)),
           Secondary = measures(fitj(oth, covs, TRUE)))

rowsA <- function(lab, f) tibble(Measure = lab,
  `Crude (95% CI)` = f("Crude"), `Primary (95% CI)` = f("Primary"), `Secondary (95% CI)` = f("Secondary"))
tabA <- bind_rows(
  tibble(Measure = "Boston (290 cases / 9,660)"),
  rowsA("  Both vs neither, aRR", function(m) one_ci(mb[[m]]$`Both MVM and AI`)),
  rowsA("  Ratio of RRs", function(m) one_ci(mb[[m]]$ratio)),
  rowsA("  RERI", function(m) one_ci(mb[[m]]$reri, log = FALSE)),
  tibble(Measure = "Other 11 sites (410 cases / 31,040)"),
  rowsA("  Both vs neither, aRR", function(m) one_ci(mo[[m]]$`Both MVM and AI`)),
  rowsA("  Ratio of RRs", function(m) one_ci(mo[[m]]$ratio)),
  rowsA("  RERI", function(m) one_ci(mo[[m]]$reri, log = FALSE)),
  tibble(Measure = "Heterogeneity, Boston vs other 11 sites"),
  rowsA("  Ratio of aRRs, MVM only", function(m) diff_ci(mb[[m]]$`MVM only`, mo[[m]]$`MVM only`)),
  rowsA("  Ratio of aRRs, AI only", function(m) diff_ci(mb[[m]]$`AI only`, mo[[m]]$`AI only`)),
  rowsA("  Ratio of aRRs, both vs neither", function(m) diff_ci(mb[[m]]$`Both MVM and AI`, mo[[m]]$`Both MVM and AI`)),
  rowsA("  Ratio of the two ratios of RRs (1 = same)", function(m) diff_ci(mb[[m]]$ratio, mo[[m]]$ratio)),
  rowsA("  Difference in RERI (0 = same)", function(m) diff_ci(mb[[m]]$reri, mo[[m]]$reri, log = FALSE))) %>%
  mutate(across(everything(), ~ replace_na(.x, "")), n_cases = 700, n_obs = 40700, n_imps = M)
write_csv(tabA, file.path(out_dir, "c_boston_heterogeneity_stratified.csv"))

# ------------------------------------------------------------------------------
# B. SINGLE MODEL WITH JOINT x BOSTON PRODUCT TERMS
# ------------------------------------------------------------------------------
cat("\n== B. Single model\n")
jb <- c("jb_mvm", "jb_ai", "jb_both")
specB <- list(Crude = list(cv = "boston", w = FALSE),
              Primary = list(cv = c(covs, "site"), w = TRUE),
              Secondary = list(cv = c(covs, "boston"), w = TRUE))
resB <- imap(specB, function(s, lab) {
  cat("   ", lab, "\n")
  fit_pool(paste(c("joint", jb), collapse = " + "), s$cv, s$w, keep = "^jb_")
})
gB <- function(res, g) {
  b <- res$est$b[match(jb, res$est$term)]; V <- res$V[jb, jb]
  one_ci(list(est = sum(g * b), var = as.numeric(t(g) %*% V %*% g)))
}
rowsB <- function(lab, g) tibble(Measure = lab, `Crude (95% CI)` = gB(resB$Crude, g),
  `Primary (95% CI)` = gB(resB$Primary, g), `Secondary (95% CI)` = gB(resB$Secondary, g))
tabB <- bind_rows(
  tibble(Measure = "Heterogeneity, Boston vs other 11 sites (single model, shared covariate effects)"),
  rowsB("  Ratio of aRRs, MVM only", c(1, 0, 0)),
  rowsB("  Ratio of aRRs, AI only", c(0, 1, 0)),
  rowsB("  Ratio of aRRs, both vs neither", c(0, 0, 1)),
  rowsB("  Ratio of the two ratios of RRs (1 = same)", c(-1, -1, 1))) %>%
  mutate(across(everything(), ~ replace_na(.x, "")), n_cases = 700, n_obs = 40700, n_imps = M)
write_csv(tabB, file.path(out_dir, "c_boston_heterogeneity_model.csv"))

cat("\n==== A. Stratified ====\n"); print(tabA, n = Inf, width = Inf)
cat("\n==== B. Single model ====\n"); print(tabB, n = Inf, width = Inf)
cat("\nDone\n")
