# ==============================================================================
# 62-joint-reruns-clean.R
# COPY of 52-joint-reruns.R for the rerun after 58 (implausible values set to missing).
# Only change: file names (_m50y -> _m50c, analytic_sample_corrected ->
# analytic_sample_clean, and c_ prefixes on CSV outputs), so the earlier
# results are kept for comparison. Code otherwise identical.
# ==============================================================================

# ==============================================================================
# 52-joint-reruns.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Three reruns for the joint exposure (neither, MVM only, AI only, both),
# using the current pipeline (49 -> 50 -> 51: 50 imputations, outcome in the
# imputation model, stabilised IPCW). Data prep, fit_pool() and the
# interaction code are copied unchanged from 51. The primary analysis is
# not changed.
#
#   0. Sanity checks: joint (crude, primary, secondary), any MVM and any AI
#      (primary, 3 decimals)
#   1. Suspect IH counted as IH (ih2, 930 cases), joint exposure, crude,
#      primary, secondary, with RERI, AP and ratio of RRs
#   2. Boston stratification, joint exposure. Boston: no site term.
#      Other 11 sites: site term (11 levels). Cases by joint category first,
#      flagging any category with fewer than 10 cases
#   3. RERI, AP and ratio of RRs for the main analysis, all three models
#
# Input:  data/weighted_long_m50c.RDS, data/path2all.sas7bdat
# Output: outputs/rerun0_sanity.csv, outputs/rerun1_ih2_joint.csv,
#         outputs/rerun1_ih2_interaction.csv, outputs/rerun2_boston_counts.csv,
#         outputs/rerun2_boston_joint.csv, outputs/rerun2_boston_interaction.csv,
#         outputs/rerun3_main_interaction.csv
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

# Interaction measures (copied from 51; delta method on the pooled covariance)
ci <- function(est, g, V, log = FALSE) {
  se <- sqrt(as.numeric(t(g) %*% V %*% g))
  v <- c(est, est - 1.96 * se, est + 1.96 * se); if (log) v <- exp(v)
  sprintf("%.2f (%.2f, %.2f)", v[1], v[2], v[3])
}
jt <- c("jointMVM only", "jointAI only", "jointBoth MVM and AI")
inter <- function(res) {
  b <- res$est$b[match(jt, res$est$term)]; V <- res$V[jt, jt]
  r10 <- exp(b[1]); r01 <- exp(b[2]); r11 <- exp(b[3]); reri <- r11 - r10 - r01 + 1
  c(`RERI (additive, 0 = none)` = ci(reri, c(-r10, -r01, r11), V),
    `Attributable proportion, AP (additive, 0 = none)` =
      ci(reri / r11, c(-r10 / r11, -r01 / r11, (r10 + r01 - 1) / r11), V),
    `Ratio of RRs (multiplicative, 1 = none)` = ci(b[3] - b[1] - b[2], c(-1, -1, 1), V, TRUE))
}

get <- function(res, term) if (is.null(res)) "" else res$est %>% filter(.data$term == !!term) %>% pull(rr_ci)

# Table in the requested layout
joint_table <- function(res, d, outcome, n_cases) {
  map_dfr(levels(d$joint), function(l) {
    k <- d[[outcome]][d$joint == l]
    tibble(`Exposure category` = if (l == "Neither") "Neither (reference)" else l,
           `IH % (cases/N)` = pct(sum(k), length(k)),
           `Crude RR (95% CI)`     = if (l == "Neither") "1.00 (Reference)" else get(res$Crude, paste0("joint", l)),
           `Primary aRR (95% CI)`  = if (l == "Neither") "1.00 (Reference)" else get(res$Primary, paste0("joint", l)),
           `Secondary aRR (95% CI)` = if (l == "Neither") "1.00 (Reference)" else get(res$Secondary, paste0("joint", l)),
           cases_flag = if (sum(k) < 10) "UNSTABLE: fewer than 10 cases" else "")
  })
}
inter_table <- function(res, n_cases, n_obs) {
  x <- map(res, ~ if (is.null(.x)) NULL else inter(.x))
  tibble(Measure = names(x$Primary),
         `Crude (95% CI)` = if (is.null(x$Crude)) "" else unname(x$Crude),
         `Primary (95% CI)` = unname(x$Primary),
         `Secondary (95% CI)` = if (is.null(x$Secondary)) "" else unname(x$Secondary),
         n_cases = n_cases, n_obs = n_obs)
}
run3 <- function(outcome, data = imp, specs = spec) {
  imap(specs, function(s, lab) { cat("   ", outcome, "/", lab, "\n"); fit_pool("joint", s$covs, s$w, outcome, data, "^joint") })
}

dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 0. SANITY CHECKS
# ------------------------------------------------------------------------------
cat("\n== 0. Sanity checks\n")
r_main <- run3("ih1")
r_mvm <- fit_pool("mvm2n", c(covs, "site"), TRUE, keep = "^mvm2n$")
r_ai  <- fit_pool("ain",   c(covs, "site"), TRUE, keep = "^ain$")
f3 <- function(e) sprintf("%.3f (%.3f, %.3f)", e$rr, e$lci, e$uci)
sanity <- bind_rows(
  imap_dfr(r_main, function(r, lab) r$est %>% transmute(model = lab, term = sub("^joint", "", term),
                                                         rr_ci_2dp = rr_ci, rr_ci_3dp = f3(r$est))),
  tibble(model = "Primary", term = "Any MVM vs none", rr_ci_2dp = r_mvm$est$rr_ci, rr_ci_3dp = f3(r_mvm$est)),
  tibble(model = "Primary", term = "Any AI vs none",  rr_ci_2dp = r_ai$est$rr_ci,  rr_ci_3dp = f3(r_ai$est))) %>%
  mutate(n_cases = 700, n_obs = 40700, n_imps = M)
write_csv(sanity, file.path(out_dir, "c_rerun0_sanity.csv"))
print(sanity, n = Inf, width = Inf)

# ------------------------------------------------------------------------------
# 1. SUSPECT IH AS IH (ih2)
# ------------------------------------------------------------------------------
cat("\n== 1. Suspect IH (ih2)\n")
r_ih2 <- run3("ih2")
t1 <- joint_table(r_ih2, d1, "ih2") %>% mutate(n_cases = sum(d1$ih2), n_obs = nrow(d1))
i1 <- inter_table(r_ih2, sum(d1$ih2), nrow(d1))
write_csv(t1, file.path(out_dir, "c_rerun1_ih2_joint.csv"))
write_csv(i1, file.path(out_dir, "c_rerun1_ih2_interaction.csv"))

# ------------------------------------------------------------------------------
# 2. BOSTON STRATIFICATION
# ------------------------------------------------------------------------------
cat("\n== 2. Boston stratification\n")
bos <- lapply(imp, function(d) d[d$site == "5", ])
oth <- lapply(imp, function(d) { d <- d[d$site != "5", ]; d$site <- droplevels(d$site); d })
stopifnot(sum(bos[[1]]$ih1) == 290, sum(oth[[1]]$ih1) == 410, nlevels(oth[[1]]$site) == 11)

counts <- bind_rows(Boston = bos[[1]], `Other 11 sites` = oth[[1]], .id = "stratum") %>%
  group_by(stratum, joint) %>%
  summarise(cases = sum(ih1), n = n(), .groups = "drop") %>%
  mutate(`IH % (cases/N)` = pct(cases, n),
         flag = if_else(cases < 10, "UNSTABLE: fewer than 10 cases", ""))
write_csv(counts, file.path(out_dir, "c_rerun2_boston_counts.csv"))
print(counts, width = Inf)

# Boston: no site term in any model, so primary and secondary are the same model
spec_bos <- list(Crude = spec$Crude, Primary = list(covs = covs, w = TRUE))
r_bos <- run3("ih1", bos, spec_bos); r_bos$Secondary <- NULL
r_oth <- run3("ih1", oth)
t2 <- bind_rows(
  joint_table(r_bos, bos[[1]], "ih1") %>% mutate(stratum = "Boston (no site term)",
    `Secondary aRR (95% CI)` = if_else(`Exposure category` == "Neither (reference)", "1.00 (Reference)", "Same as primary (no site term)"),
    n_cases = sum(bos[[1]]$ih1), n_obs = nrow(bos[[1]])),
  joint_table(r_oth, oth[[1]], "ih1") %>% mutate(stratum = "Other 11 sites (primary with 11-level site term)",
    n_cases = sum(oth[[1]]$ih1), n_obs = nrow(oth[[1]]))) %>%
  select(stratum, everything())
i2 <- bind_rows(
  inter_table(r_bos, sum(bos[[1]]$ih1), nrow(bos[[1]])) %>% mutate(stratum = "Boston (no site term)",
    `Secondary (95% CI)` = "Same as primary (no site term)"),
  inter_table(r_oth, sum(oth[[1]]$ih1), nrow(oth[[1]])) %>% mutate(stratum = "Other 11 sites")) %>%
  select(stratum, everything())
write_csv(t2, file.path(out_dir, "c_rerun2_boston_joint.csv"))
write_csv(i2, file.path(out_dir, "c_rerun2_boston_interaction.csv"))

# ------------------------------------------------------------------------------
# 3. INTERACTION, MAIN ANALYSIS, ALL THREE MODELS
# ------------------------------------------------------------------------------
t3 <- joint_table(r_main, d1, "ih1") %>% mutate(n_cases = 700, n_obs = 40700)
i3 <- inter_table(r_main, 700, 40700)
write_csv(t3, file.path(out_dir, "c_rerun3_main_joint.csv"))
write_csv(i3, file.path(out_dir, "c_rerun3_main_interaction.csv"))

cat("\n==== RERUN 1 (ih2) ====\n"); print(t1, width = Inf); print(i1, width = Inf)
cat("\n==== RERUN 2 (Boston) ====\n"); print(t2, width = Inf); print(i2, width = Inf)
cat("\n==== RERUN 3 (main) ====\n"); print(t3, width = Inf); print(i3, width = Inf)
cat("\nDone\n")
