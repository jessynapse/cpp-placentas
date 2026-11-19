# ==============================================================================
# 74-ga-implausible-sensitivity.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# David Margolis: "What happens if you drop the few that were improbable?"
# Gestational age is not in the imputation model, the censoring weights or any
# outcome model. This refits the 9-covariate primary model (71 to 73) after
# excluding examined infants whose gestational age is improbable:
#   outside 20 to 44 completed weeks, or implausible for birthweight
#   (under 28 weeks with birthweight over 2,500 g, or 37 weeks or more with
#   birthweight under 1,000 g, birthweight under 300 g or over 6,000 g ignored)
# Same 50 imputations, same 9-covariate censoring weights (ipw9), same GEE.
# Two exclusion sets: (1) improbable gestational age, (2) improbable plus
# missing gestational age.
#
# Input:  data/weighted_long_m50y.RDS
# Output: outputs/ga_implausible_sensitivity.csv, outputs/ga_implausible_counts.csv
# Cache:  data/ga_sensitivity_fits_74.RDS (git-ignored)
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
test_n  <- as.integer(Sys.getenv("R74_TEST", "0"))
if (test_n > 0) out_dir <- tempdir()
n_cores <- min(4, detectCores())
f2 <- function(r, l, u) sprintf("%.2f (%.2f, %.2f)", r, l, u)

# ------------------------------------------------------------------------------
# DATA (as 71 and 72)
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
  joint   = factor(case_when(mvm2n == 0 & ain == 0 ~ "Neither", mvm2n == 1 & ain == 0 ~ "MVM only",
                             mvm2n == 0 & ain == 1 ~ "AI only", TRUE ~ "Both MVM and AI"),
                   levels = c("Neither", "MVM only", "AI only", "Both MVM and AI"))) %>%
  mutate(across(c(dm, chronic_htn, infant_sex), ~ factor(as.character(.x))))

red9 <- c("age", "race", "educ", "income", "marital", "smoking", "parity", "infant_sex", "site")
lab  <- c(age = "Maternal age", race = "Race and ethnicity", educ = "Education", income = "Annual family income",
          marital = "Marital status", smoking = "Smoking", parity = "Parity", infant_sex = "Infant sex",
          site = "Study site")

all_imp <- w %>% filter(.imp > 0) %>% prep() %>%
  select(.imp, MOMID, ih1, exclusion, ih_obs, site, age, race, educ, income, marital, smoking, parity,
         infant_sex, mvm2, ai, joint, gest_age, birthweight) %>%
  split(.$.imp)
rm(w); invisible(gc())
M <- length(all_imp)
stopifnot(test_n > 0 || M == 50, all(sapply(all_imp, nrow) == 45268))

rw <- function(d, covs) {    # stabilised IPCW, as 50, 71 and 72
  fd <- as.formula(paste("ih_obs ~ mvm2 + ai +", paste(covs, collapse = " + ")))
  pd <- predict(glm(fd, data = d, family = binomial(link = "logit"), na.action = na.exclude), newdata = d, type = "response")
  pn <- predict(glm(ih_obs ~ 1, data = d, family = binomial(link = "logit")), newdata = d, type = "response")
  as.numeric(pn / pd)
}
imp <- lapply(all_imp, function(d) { d$ipw9 <- rw(d, red9); d %>% filter(!is.na(ih1), exclusion == 0) %>% arrange(MOMID) })
rm(all_imp); invisible(gc())
d1 <- imp[[1]]
stopifnot(nrow(d1) == 40700, sum(d1$ih1) == 700)
if (test_n == 0) {
  w71 <- read_csv(file.path(out_dir, "reduced9_weights.csv"), show_col_types = FALSE)
  stopifnot(abs(max(d1$ipw9) - w71$max[w71$weights == "ipw9"]) < 1e-8)   # same weights as 71 and 72
}

# ------------------------------------------------------------------------------
# FIT AND POOL (as 71), all models weighted with ipw9
# ------------------------------------------------------------------------------
cache_f <- file.path(in_dir, "ga_sensitivity_fits_74.RDS")
cache <- if (test_n == 0 && file.exists(cache_f)) readRDS(cache_f) else list()
jt <- c(`Maternal vascular malperfusion only` = "jointMVM only", `Acute inflammation only` = "jointAI only",
        Both = "jointBoth MVM and AI")
fit_pool <- function(key, covs, data = imp) {
  if (!is.null(cache[[key]])) return(cache[[key]])
  cat("Fitting", key, "\n")
  f <- as.formula(paste("ih1 ~", paste(c("joint", covs), collapse = " + ")))
  one <- function(d) tryCatch({
    m <- geeglm(f, data = d, family = poisson(link = "log"), id = MOMID, corstr = "exchangeable", weights = ipw9)
    list(b = coef(m), V = vcov(m)) }, error = function(e) NULL)
  fits <- mclapply(data, one, mc.cores = n_cores)
  for (i in which(sapply(fits, is.null))) fits[[i]] <- one(data[[i]])
  fits <- Filter(Negate(is.null), fits)
  stopifnot(length(fits) == length(data))
  m <- length(fits)
  B <- do.call(rbind, lapply(fits, `[[`, "b"))[, jt, drop = FALSE]
  W <- Reduce(`+`, lapply(fits, function(x) x$V[jt, jt, drop = FALSE])) / m
  V <- W + (1 + 1/m) * cov(B)
  b <- colMeans(B); se <- sqrt(diag(V))
  res <- tibble(exposure = names(jt), rr = exp(b), lci = exp(b - 1.96 * se), uci = exp(b + 1.96 * se))
  cache[[key]] <<- res
  if (test_n == 0) saveRDS(cache, cache_f)
  res
}


# ------------------------------------------------------------------------------
# IMPROBABLE GESTATIONAL AGE (exposure, outcome and GA are never imputed)
# ------------------------------------------------------------------------------
flag <- function(d) {
  bw <- ifelse(d$birthweight < 300 | d$birthweight > 6000, NA, d$birthweight)
  g <- d$gest_age
  !is.na(g) & (g < 20 | g > 44 | (!is.na(bw) & ((g < 28 & bw > 2500) | (g >= 37 & bw < 1000))))
}
imp <- lapply(imp, function(d) { d$ga_bad <- flag(d); d })
d1 <- imp[[1]]
stopifnot(all(sapply(imp, function(d) identical(d$ga_bad, d1$ga_bad))))
pc <- function(k, n) sprintf("%.1f%% (%s/%s)", 100 * k / n, format(k, big.mark = ","), format(n, big.mark = ","))
counts <- tibble(
  item = c("Gestational age under 20 weeks", "Gestational age over 44 weeks",
           "Implausible for birthweight (within 20 to 44 weeks)", "Improbable gestational age, total",
           "IH cases among infants with improbable gestational age", "Gestational age missing"),
  value = c(pc(sum(d1$gest_age < 20, na.rm = TRUE), 40700), pc(sum(d1$gest_age > 44, na.rm = TRUE), 40700),
            pc(sum(d1$ga_bad & d1$gest_age >= 20 & d1$gest_age <= 44), 40700), pc(sum(d1$ga_bad), 40700),
            pc(sum(d1$ih1[d1$ga_bad]), 700), pc(sum(is.na(d1$gest_age)), 40700)))
print(counts)
by_cat <- d1 %>% group_by(joint) %>% summarise(n = n(), improbable = sum(ga_bad), cases_dropped = sum(ih1[ga_bad]), .groups = "drop")
print(by_cat)

keep1 <- lapply(imp, function(d) d[!d$ga_bad, ])
keep2 <- lapply(imp, function(d) d[!d$ga_bad & !is.na(d$gest_age), ])

full <- fit_pool("full", red9)
both_full <- with(full[full$exposure == "Both", ], f2(rr, lci, uci))
cat("Full 9-covariate model, both patterns:", both_full, "(target 1.36 (1.04, 1.77))\n")
if (test_n == 0 && both_full != "1.36 (1.04, 1.77)") stop("Primary model not reproduced. Stopping.")

mods <- list(
  list(name = "Primary, all examined infants", data = imp, covs = red9, key = "full"),
  list(name = "Primary, excluding improbable gestational age", data = keep1, covs = red9, key = "drop_bad"),
  list(name = "Primary, excluding improbable or missing gestational age", data = keep2, covs = red9, key = "drop_bad_missing"),
  list(name = "Crude (weighted), all examined infants", data = imp, covs = character(0), key = "crude"),
  list(name = "Crude (weighted), excluding improbable gestational age", data = keep1, covs = character(0), key = "crude_drop_bad"))
res <- map_dfr(mods, function(m) {
  r <- fit_pool(m$key, m$covs, data = m$data)
  r %>% mutate(model = m$name, n_obs = nrow(m$data[[1]]), n_cases = sum(m$data[[1]]$ih1))
}) %>%
  group_by(base = if_else(grepl("^Crude", model), "crude", "adjusted"), exposure) %>%
  mutate(pct_change_vs_all = 100 * (rr - rr[grepl("all examined", model)]) / rr[grepl("all examined", model)]) %>%
  ungroup() %>%
  transmute(model, exposure, `RR (95% CI)` = f2(rr, lci, uci), RR = round(rr, 3), lower_ci = round(lci, 3),
            upper_ci = round(uci, 3), pct_change_vs_all = round(pct_change_vs_all, 1), n_cases, n_obs, n_imps = M)
write_csv(res, file.path(out_dir, "ga_implausible_sensitivity.csv"))
write_csv(bind_rows(counts, by_cat %>% transmute(item = paste("Improbable gestational age, category:", joint),
  value = sprintf("%s of %s infants, %s IH cases", format(improbable, big.mark = ","), format(n, big.mark = ","), cases_dropped))),
  file.path(out_dir, "ga_implausible_counts.csv"))
print(res, n = Inf, width = Inf)
cat("\nDone\n")
