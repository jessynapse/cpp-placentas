# ==============================================================================
# 72-v52-tables-reduced9-primary.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Final table values for manuscript v5.2 with the 9-covariate model (71) as
# the primary analysis, in the v5.1 layout:
#   Table 2   four-category pattern, unadjusted and adjusted, RERI, AP, ratio of RRs
#   eTable 2  each pattern any vs none
#   eTable 3  MVM grade and AI stage
#   eTable 4  Boston vs other 11 sites
#   eTable 5  definite vs definite or suspect IH
#   eTable 6  E-values and maternal history of angioma bounds
#   eTable 7  12-covariate vs 9-covariate model, four-category pattern
# Plus missing covariate data for the 9 covariates (and the 12, as a check).
#
# Pooled fits are read from the 71 cache (no refitting of 71 models). The
# unadjusted column is unweighted with no covariates, as in 51 and 71, so it
# does not depend on the censoring weights. As a check, weighted unadjusted
# models are fitted here with the 12-covariate and the 9-covariate weights.
#
# Input:  data/weighted_long_m50y.RDS, data/reduced9_fits_71.RDS,
#         outputs/reduced9_comparison.csv
# Output: outputs/v5.2_table2.csv, outputs/v5.2_etable2.csv ... v5.2_etable7.csv,
#         outputs/v5.2_missing_covariates.csv, outputs/v5.2_unadjusted_weight_check.csv
# Cache:  data/v52_fits_72.RDS (git-ignored, weighted unadjusted fits)
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
n_cores <- min(4, detectCores())
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)
pct   <- function(k, n) sprintf("%.1f%% (%s/%s)", 100 * k / n, fmt_n(k), fmt_n(n))
f2    <- function(r, l, u) sprintf("%.2f (%.2f, %.2f)", r, l, u)
fx    <- function(x) f2(x[[1]], x[[2]], x[[3]])
REF   <- "1.00 (Reference)"

# ------------------------------------------------------------------------------
# DATA (as 71)
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
  mvm3    = factor(as.character(mvm3), levels = c("None", "Low", "High")),
  ai3     = factor(as.character(ai3),  levels = c("None", "Low", "High")),
  joint   = factor(case_when(mvm2n == 0 & ain == 0 ~ "Neither", mvm2n == 1 & ain == 0 ~ "MVM only",
                             mvm2n == 0 & ain == 1 ~ "AI only", TRUE ~ "Both MVM and AI"),
                   levels = c("Neither", "MVM only", "AI only", "Both MVM and AI"))) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), ~ factor(as.character(.x))))

full12 <- c("age", "bmi", "race", "educ", "income", "marital", "smoking", "parity",
            "dm", "chronic_htn", "infant_sex", "site")
red9   <- setdiff(full12, c("bmi", "dm", "chronic_htn"))

# Missing covariate data, observed data (.imp == 0), 40,700 examined infants
obs0 <- w %>% filter(.imp == 0, !is.na(ih1), exclusion == 0)
stopifnot(nrow(obs0) == 40700)
lab <- c(age = "Maternal age", bmi = "Pre-pregnancy BMI", race = "Race and ethnicity", educ = "Education",
         income = "Annual family income", marital = "Marital status", smoking = "Smoking", parity = "Parity",
         dm = "Pre-pregnancy diabetes", chronic_htn = "Chronic hypertension", infant_sex = "Infant sex",
         site = "Study site")
any_miss <- function(v) rowSums(is.na(as.data.frame(obs0)[, v])) > 0
miss <- bind_rows(
  tibble(set = "9 covariates", covariate = lab[red9], missing = map_chr(red9, ~ pct(sum(is.na(obs0[[.x]])), 40700))),
  tibble(set = "9 covariates", covariate = "At least one of the 9", missing = pct(sum(any_miss(red9)), 40700)),
  tibble(set = "12 covariates", covariate = lab[full12], missing = map_chr(full12, ~ pct(sum(is.na(obs0[[.x]])), 40700))),
  tibble(set = "12 covariates", covariate = "At least one of the 12", missing = pct(sum(any_miss(full12)), 40700)))
print(miss, n = Inf)
stopifnot(sprintf("%.1f", 100 * mean(any_miss(full12))) == "14.6")

all_imp <- w %>% filter(.imp > 0) %>% prep() %>%
  select(.imp, MOMID, ih1, ih2, exclusion, ih_obs, ipw, site, age, bmi, race, educ, income, marital,
         smoking, parity, dm, chronic_htn, infant_sex, mvm2, ai, mvm2n, ain, mvm3, ai3, joint) %>%
  split(.$.imp)
rm(w, obs0); invisible(gc())
stopifnot(length(all_imp) == 50)

rw <- function(d, covs) {    # stabilised IPCW, as 50 and 71
  fd <- as.formula(paste("ih_obs ~ mvm2 + ai +", paste(covs, collapse = " + ")))
  pd <- predict(glm(fd, data = d, family = binomial(link = "logit"), na.action = na.exclude), newdata = d, type = "response")
  pn <- predict(glm(ih_obs ~ 1, data = d, family = binomial(link = "logit")), newdata = d, type = "response")
  as.numeric(pn / pd)
}
imp <- lapply(all_imp, function(d) { d$ipw9 <- rw(d, red9); d %>% filter(!is.na(ih1), exclusion == 0) %>% arrange(MOMID) })
rm(all_imp); invisible(gc())
d1 <- imp[[1]]
stopifnot(nrow(d1) == 40700, sum(d1$ih1) == 700, sum(d1$ih2) == 930)
w71 <- read_csv(file.path(out_dir, "reduced9_weights.csv"), show_col_types = FALSE)
stopifnot(abs(max(d1$ipw9) - w71$max[w71$weights == "ipw9"]) < 1e-8)   # same weights as 71

# ------------------------------------------------------------------------------
# POOLED FITS: 71 cache (read only) plus weighted unadjusted fits
# ------------------------------------------------------------------------------
c71 <- readRDS(file.path(in_dir, "reduced9_fits_71.RDS"))
cache_f <- file.path(in_dir, "v52_fits_72.RDS")
c72 <- if (file.exists(cache_f)) readRDS(cache_f) else list()
fit_pool <- function(key, expo, wcol, data = imp) {
  if (!is.null(c72[[key]])) return(c72[[key]])
  cat("Fitting", key, "\n")
  f <- as.formula(paste("ih1 ~", expo))
  one <- function(d) { d$.w <- d[[wcol]]
    m <- geeglm(f, data = d, family = poisson(link = "log"), id = MOMID, corstr = "exchangeable", weights = .w)
    list(b = coef(m), V = vcov(m)) }
  fits <- mclapply(data, one, mc.cores = n_cores)
  stopifnot(!any(sapply(fits, inherits, "try-error")))
  m <- length(fits)
  B <- do.call(rbind, lapply(fits, `[[`, "b")); terms <- grep(paste0("^", expo), colnames(B), value = TRUE)
  B <- B[, terms, drop = FALSE]
  W <- Reduce(`+`, lapply(fits, function(x) x$V[terms, terms, drop = FALSE])) / m
  res <- list(b = colMeans(B), V = W + (1 + 1/m) * cov(B), n_imps = m)
  c72[[key]] <<- res; saveRDS(c72, cache_f); res
}
lin <- function(fit, g) { e <- sum(g * fit$b); s <- sqrt(as.numeric(t(g) %*% fit$V %*% g))
  c(rr = exp(e), lci = exp(e - 1.96 * s), uci = exp(e + 1.96 * s)) }
term_rr <- function(fit, term) lin(fit, as.numeric(names(fit$b) == term))
inter <- function(fit) {
  r <- exp(fit$b[c("jointMVM only", "jointAI only", "jointBoth MVM and AI")]); V <- fit$V[names(r), names(r)]
  dm <- function(est, g) { s <- sqrt(as.numeric(t(g) %*% V %*% g)); c(est, est - 1.96 * s, est + 1.96 * s) }
  reri <- r[3] - r[1] - r[2] + 1
  g <- c(-1, -1, 1); e <- sum(g * log(r)); s <- sqrt(as.numeric(t(g) %*% V %*% g))
  list(RERI = dm(reri, c(-r[1], -r[2], r[3])),
       AP = dm(reri / r[3], c(-r[1] / r[3], -r[2] / r[3], (r[1] + r[2] - 1) / r[3])),
       `Ratio of RRs` = c(exp(e), exp(e - 1.96 * s), exp(e + 1.96 * s)))
}
ev    <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
ev_ci <- function(l, u) ifelse(l <= 1 & u >= 1, 1, ifelse(l > 1, ev(l), ev(u)))
bound <- function(rr, ud = 1.64) rr * (ud - 1) / (ud - rr)

jt <- c(`Maternal vascular malperfusion only` = "jointMVM only", `Acute inflammation only` = "jointAI only",
        Both = "jointBoth MVM and AI")
jl <- c(Neither = "Neither", `Maternal vascular malperfusion only` = "MVM only",
        `Acute inflammation only` = "AI only", Both = "Both MVM and AI")
cnt <- function(x, y = d1$ih1) pct(sum(y[x]), sum(x))
P <- c71[c("9_joint", "9_mvm2n", "9_ain", "9_mvm3", "9_ai3", "9_bos", "9_oth", "9_ih2")]
names(P) <- sub("^9_", "", names(P))
U <- c71[c("crude_joint", "crude_mvm2n", "crude_ain", "crude_mvm3", "crude_ai3")]
names(U) <- sub("^crude_", "", names(U))

# Check against 71's comparison file
cmp71 <- read_csv(file.path(out_dir, "reduced9_comparison.csv"), show_col_types = FALSE)
stopifnot(fx(term_rr(P$joint, jt[["Both"]])) == "1.36 (1.04, 1.77)",
          fx(term_rr(P$joint, jt[["Both"]])) == cmp71$`Comparison model`[cmp71$table == "Table 2" & cmp71$exposure == "Both" & grepl("^9", cmp71$model)])

# ------------------------------------------------------------------------------
# CHECK: does a weighted unadjusted model change with the 9-covariate weights?
# ------------------------------------------------------------------------------
wchk <- map_dfr(c(joint = "joint", mvm2n = "mvm2n", ain = "ain", mvm3 = "mvm3", ai3 = "ai3"), function(e) {
  f12 <- fit_pool(paste0("w12_", e), e, "ipw"); f9 <- fit_pool(paste0("w9_", e), e, "ipw9")
  tibble(term = names(U[[e]]$b),
         `Unweighted (manuscript column)` = map_chr(term, ~ fx(term_rr(U[[e]], .x))),
         `Weighted, 12-covariate weights` = map_chr(term, ~ fx(term_rr(f12, .x))),
         `Weighted, 9-covariate weights` = map_chr(term, ~ fx(term_rr(f9, .x))))
})
wi <- map_dfr(list(Unweighted = U$joint, w12 = c72$w12_joint, w9 = c72$w9_joint), ~ map_chr(inter(.x), fx), .id = "model")
wchk <- bind_rows(wchk, tibble(term = c("RERI", "AP", "Ratio of RRs"),
  `Unweighted (manuscript column)` = unlist(wi[1, -1]), `Weighted, 12-covariate weights` = unlist(wi[2, -1]),
  `Weighted, 9-covariate weights` = unlist(wi[3, -1])))
print(wchk, n = Inf, width = Inf)

# ------------------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------------------
t2 <- bind_rows(
  tibble(Exposure = names(jl), `IH, % (cases/N)` = map_chr(jl, ~ cnt(d1$joint == .x)),
         `Unadjusted RR (95% CI)` = c(REF, map_chr(jt, ~ fx(term_rr(U$joint, .x)))),
         `aRR (95% CI)` = c(REF, map_chr(jt, ~ fx(term_rr(P$joint, .x))))),
  tibble(Exposure = c("RERI (additive, 0 = none)", "AP (additive, 0 = none)", "Ratio of RRs (multiplicative, 1 = none)"),
         `IH, % (cases/N)` = "", `Unadjusted RR (95% CI)` = map_chr(inter(U$joint), fx),
         `aRR (95% CI)` = map_chr(inter(P$joint), fx)))

et2 <- bind_rows(
  tibble(Block = "Maternal vascular malperfusion, any vs none", Level = c("None", "Any"),
         `IH, % (cases/N)` = c(cnt(d1$mvm2n == 0), cnt(d1$mvm2n == 1)),
         `Unadjusted RR (95% CI)` = c(REF, fx(term_rr(U$mvm2n, "mvm2n"))), `aRR (95% CI)` = c(REF, fx(term_rr(P$mvm2n, "mvm2n")))),
  tibble(Block = "Acute inflammation, any vs none", Level = c("None", "Any"),
         `IH, % (cases/N)` = c(cnt(d1$ain == 0), cnt(d1$ain == 1)),
         `Unadjusted RR (95% CI)` = c(REF, fx(term_rr(U$ain, "ain"))), `aRR (95% CI)` = c(REF, fx(term_rr(P$ain, "ain")))))

et3 <- bind_rows(
  tibble(Block = "Maternal vascular malperfusion grade", Level = c("None", "Low (score 1 to 2)", "High (score 3 or more)"),
         `IH, % (cases/N)` = map_chr(c("None", "Low", "High"), ~ cnt(d1$mvm3 == .x)),
         `Unadjusted RR (95% CI)` = c(REF, fx(term_rr(U$mvm3, "mvm3Low")), fx(term_rr(U$mvm3, "mvm3High"))),
         `aRR (95% CI)` = c(REF, fx(term_rr(P$mvm3, "mvm3Low")), fx(term_rr(P$mvm3, "mvm3High")))),
  tibble(Block = "Acute inflammation stage", Level = c("None", "Low", "High"),
         `IH, % (cases/N)` = map_chr(c("None", "Low", "High"), ~ cnt(d1$ai3 == .x)),
         `Unadjusted RR (95% CI)` = c(REF, fx(term_rr(U$ai3, "ai3Low")), fx(term_rr(U$ai3, "ai3High"))),
         `aRR (95% CI)` = c(REF, fx(term_rr(P$ai3, "ai3Low")), fx(term_rr(P$ai3, "ai3High")))))

bd <- d1[d1$site == "5", ]; od <- d1[d1$site != "5", ]
sc <- function(dd, x, y = "ih1") pct(sum(dd[[y]][dd$joint == x]), sum(dd$joint == x))
et4 <- tibble(Exposure = c("All infants", names(jl)),
  `Boston, IH % (cases/N)` = c(pct(sum(bd$ih1), nrow(bd)), map_chr(jl, ~ sc(bd, .x))),
  `Boston, aRR (95% CI)` = c("", REF, map_chr(jt, ~ fx(term_rr(P$bos, .x)))),
  `Other 11 sites, IH % (cases/N)` = c(pct(sum(od$ih1), nrow(od)), map_chr(jl, ~ sc(od, .x))),
  `Other 11 sites, aRR (95% CI)` = c("", REF, map_chr(jt, ~ fx(term_rr(P$oth, .x)))))

et5 <- tibble(Exposure = c("All infants", names(jl)),
  `Definite IH, IH % (cases/N)` = c(pct(700, 40700), map_chr(jl, ~ cnt(d1$joint == .x))),
  `Definite IH, aRR (95% CI)` = c("", REF, map_chr(jt, ~ fx(term_rr(P$joint, .x)))),
  `Definite or suspect IH, IH % (cases/N)` = c(pct(930, 40700), map_chr(jl, ~ cnt(d1$joint == .x, d1$ih2))),
  `Definite or suspect IH, aRR (95% CI)` = c("", REF, map_chr(jt, ~ fx(term_rr(P$ih2, .x)))))

et6 <- imap_dfr(jt, function(t, e) { x <- term_rr(P$joint, t)
  tibble(Exposure = paste(e, "vs neither"), `aRR (95% CI)` = fx(x),
         `E-value, point estimate` = sprintf("%.2f", ev(x[["rr"]])), `E-value, CI bound` = sprintf("%.2f", ev_ci(x[["lci"]], x[["uci"]]))) })
xb <- term_rr(P$joint, jt[["Both"]])
et6b <- tibble(Quantity = c("RR_EU needed to fully explain the aRR for both patterns",
                            "RR_EU needed to move the lower confidence limit to 1"),
               `Based on` = c(sprintf("aRR %.2f", xb[["rr"]]), sprintf("lower limit %.2f", xb[["lci"]])),
               Value = sprintf("%.2f", c(bound(xb[["rr"]]), bound(xb[["lci"]]))))

j12 <- c71[["12_joint"]]
et7 <- tibble(Exposure = c(names(jl), "RERI (additive, 0 = none)", "AP (additive, 0 = none)", "Ratio of RRs (multiplicative, 1 = none)"),
  `IH, % (cases/N)` = c(map_chr(jl, ~ cnt(d1$joint == .x)), "", "", ""),
  `12 covariates, aRR (95% CI)` = c(REF, map_chr(jt, ~ fx(term_rr(j12, .x))), map_chr(inter(j12), fx)),
  `9 covariates (primary), aRR (95% CI)` = c(REF, map_chr(jt, ~ fx(term_rr(P$joint, .x))), map_chr(inter(P$joint), fx)))
stopifnot(et7$`12 covariates, aRR (95% CI)`[4] == "1.37 (1.05, 1.79)")

tabs <- list(table2 = t2, etable2 = et2, etable3 = et3, etable4 = et4, etable5 = et5, etable6 = et6,
             etable6_angioma = et6b, etable7 = et7, missing_covariates = miss, unadjusted_weight_check = wchk)
iwalk(tabs, ~ write_csv(.x, file.path(out_dir, paste0("v5.2_", .y, ".csv"))))
iwalk(tabs, ~ { cat("\n==", .y, "\n"); print(.x, n = Inf, width = Inf) })
cat("\nDone\n")
