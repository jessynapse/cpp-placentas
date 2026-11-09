# ==============================================================================
# 64-interaction-knol-clean.R
# COPY of 54-interaction-knol-vanderweele.R for the rerun after 58 (implausible values set to missing).
# Only change: file names (_m50y -> _m50c, analytic_sample_corrected ->
# analytic_sample_clean, and c_ prefixes on CSV outputs), so the earlier
# results are kept for comparison. Code otherwise identical.
# ==============================================================================

# ==============================================================================
# 54-interaction-knol-vanderweele.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Full interaction table for the joint MVM x AI exposure, in the format of
# Knol and VanderWeele (2012, Int J Epidemiol 41:514-520):
#   - cases/N and RR (95% CI) for each of the four cells vs one reference
#   - RR of each exposure within strata of the other
#       MVM when AI absent  = exp(b_MVMonly)
#       MVM when AI present = exp(b_both - b_AIonly)
#       AI when MVM absent  = exp(b_AIonly)
#       AI when MVM present = exp(b_both - b_MVMonly)
#   - RERI, AP and ratio of RRs (as 51 and 52)
# Crude, primary and secondary models. Same data prep, weights, fit_pool()
# and delta-method code as 51 and 52 (pooled covariance by Rubin's rules).
# Primary analysis unchanged. No p-values.
#
# Input:  data/weighted_long_m50c.RDS
# Output: outputs/interaction_knol_cells.csv, outputs/interaction_knol_strata.csv,
#         outputs/interaction_knol_measures.csv
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
n_cores <- min(4, detectCores())   # set to 1 on Windows

fmt_ci <- function(rr, l, u) sprintf("%.2f (%.2f, %.2f)", rr, l, u)
fmt_n  <- function(x) format(x, big.mark = ",", trim = TRUE)
pct    <- function(k, n) sprintf("%.1f%% (%s/%s)", 100 * k / n, fmt_n(k), fmt_n(n))

# ------------------------------------------------------------------------------
# DATA (as 51)
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
                     levels = c("Neither", "MVM only", "AI only", "Both MVM and AI"))
  ) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), ~ factor(as.character(.x)))) %>%
  select(.imp, MOMID, ih1, ih2, ipw, site, age, bmi, race, educ, income,
         marital, smoking, parity, dm, chronic_htn, infant_sex, mvm2n, ain, joint) %>%
  arrange(MOMID) %>%
  split(.$.imp)
rm(w); invisible(gc())
d1 <- imp[[1]]
M <- length(imp)
stopifnot(nrow(d1) == 40700, sum(d1$ih1) == 700, sum(d1$ih2) == 930, M == 50)

# ------------------------------------------------------------------------------
# FIT AND POOL (copied from 51)
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


# ------------------------------------------------------------------------------
# FIT (joint exposure, three models)
# ------------------------------------------------------------------------------
jt <- c("jointMVM only", "jointAI only", "jointBoth MVM and AI")
res <- imap(spec, function(s, lab) { cat("   ", lab, "\n"); fit_pool("joint", s$covs, s$w, keep = "^joint") })

# Linear combination of the three log RRs, with a delta-method CI
lc <- function(r, g, log = TRUE) {
  b <- r$est$b[match(jt, r$est$term)]; V <- r$V[jt, jt]
  est <- sum(g * b); se <- sqrt(as.numeric(t(g) %*% V %*% g))
  v <- c(est, est - 1.96 * se, est + 1.96 * se); if (log) v <- exp(v)
  sprintf("%.2f (%.2f, %.2f)", v[1], v[2], v[3])
}
nonlin <- function(r) {
  b <- r$est$b[match(jt, r$est$term)]; V <- r$V[jt, jt]
  r10 <- exp(b[1]); r01 <- exp(b[2]); r11 <- exp(b[3]); reri <- r11 - r10 - r01 + 1
  f <- function(est, g) { se <- sqrt(as.numeric(t(g) %*% V %*% g)); sprintf("%.2f (%.2f, %.2f)", est, est - 1.96 * se, est + 1.96 * se) }
  c(RERI = f(reri, c(-r10, -r01, r11)),
    AP = f(reri / r11, c(-r10 / r11, -r01 / r11, (r10 + r01 - 1) / r11)))
}

# ------------------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------------------
cell <- function(m, a) { k <- d1$ih1[d1$mvm2n == m & d1$ain == a]; sprintf("%.1f%% (%s/%s)", 100 * mean(k), fmt_n(sum(k)), fmt_n(length(k))) }
g <- list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
cells <- tibble(
  `MVM` = c("None", "None", "Any", "Any"), `AI` = c("None", "Any", "None", "Any"),
  `IH % (cases/N)` = c(cell(0, 0), cell(0, 1), cell(1, 0), cell(1, 1)),
  `Crude RR (95% CI)`      = c("1.00 (Reference)", lc(res$Crude, g[[2]]),     lc(res$Crude, g[[1]]),     lc(res$Crude, g[[3]])),
  `Primary aRR (95% CI)`   = c("1.00 (Reference)", lc(res$Primary, g[[2]]),   lc(res$Primary, g[[1]]),   lc(res$Primary, g[[3]])),
  `Secondary aRR (95% CI)` = c("1.00 (Reference)", lc(res$Secondary, g[[2]]), lc(res$Secondary, g[[1]]), lc(res$Secondary, g[[3]])),
  n_cases = 700, n_obs = 40700)

strata <- tibble(
  `Effect` = c("Any MVM vs none, among placentas without AI", "Any MVM vs none, among placentas with AI",
               "Any AI vs none, among placentas without MVM", "Any AI vs none, among placentas with MVM"),
  gr = list(c(1, 0, 0), c(0, -1, 1), c(0, 1, 0), c(-1, 0, 1))) %>%
  mutate(`Crude RR (95% CI)`      = map_chr(gr, ~ lc(res$Crude, .x)),
         `Primary aRR (95% CI)`   = map_chr(gr, ~ lc(res$Primary, .x)),
         `Secondary aRR (95% CI)` = map_chr(gr, ~ lc(res$Secondary, .x)),
         n_cases = 700, n_obs = 40700) %>% select(-gr)

nl <- map(res, nonlin)
measures <- tibble(
  Measure = c("RERI (additive, 0 = none)", "Attributable proportion, AP (additive, 0 = none)",
              "Ratio of RRs (multiplicative, 1 = none)"),
  `Crude (95% CI)`     = c(nl$Crude[["RERI"]],     nl$Crude[["AP"]],     lc(res$Crude, c(-1, -1, 1))),
  `Primary (95% CI)`   = c(nl$Primary[["RERI"]],   nl$Primary[["AP"]],   lc(res$Primary, c(-1, -1, 1))),
  `Secondary (95% CI)` = c(nl$Secondary[["RERI"]], nl$Secondary[["AP"]], lc(res$Secondary, c(-1, -1, 1))),
  n_cases = 700, n_obs = 40700)

dir.create(out_dir, showWarnings = FALSE)
write_csv(cells,    file.path(out_dir, "c_interaction_knol_cells.csv"))
write_csv(strata,   file.path(out_dir, "c_interaction_knol_strata.csv"))
write_csv(measures, file.path(out_dir, "c_interaction_knol_measures.csv"))
print(cells, width = Inf); print(strata, width = Inf); print(measures, width = Inf)
cat("\nDone\n")
