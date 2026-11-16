# ==============================================================================
# 73-change-in-estimate.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Committee chair: "Is there one or two covariates that are responsible for
# the largest change in RR?" Change-in-estimate analysis with the 9-covariate
# model (71, 72) as the primary model.
#
# Pipeline as 72: 40,700 examined infants (700 definite IH), 50 imputations
# pooled with Rubin's rules, modified Poisson GEE, exchangeable, clustered on
# MOMID, stabilised IPCW from both patterns plus the 9 covariates (ipw9, as
# 71). Every model below uses these same weights, so changes reflect
# covariates only, not weighting.
#
# Models (four-category pattern, neither = reference):
#   Crude       no covariates, weighted
#   Add one     crude + each covariate alone (9 models)
#   Drop one    full minus each covariate (9 models)
#   Full        all 9 covariates (must reproduce 1.36 (1.04, 1.77) for both)
# Add-one change  = (RR_crude+cov - RR_crude) / RR_crude x 100
# Drop-one change = (RR_full-cov - RR_full) / RR_full x 100
#
# Input:  data/weighted_long_m50y.RDS
# Output: outputs/change_in_estimate.csv
# Cache:  data/change_in_estimate_fits_73.RDS (git-ignored)
# Test:   R73_TEST=3 uses 3 imputations and writes nothing to outputs
# ==============================================================================

library(tidyverse)
library(geepack)
library(parallel)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
test_n  <- as.integer(Sys.getenv("R73_TEST", "0"))
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
         infant_sex, mvm2, ai, joint) %>%
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
cache_f <- file.path(in_dir, "change_in_estimate_fits_73.RDS")
cache <- if (test_n == 0 && file.exists(cache_f)) readRDS(cache_f) else list()
jt <- c(`Maternal vascular malperfusion only` = "jointMVM only", `Acute inflammation only` = "jointAI only",
        Both = "jointBoth MVM and AI")
fit_pool <- function(key, covs) {
  if (!is.null(cache[[key]])) return(cache[[key]])
  cat("Fitting", key, "\n")
  f <- as.formula(paste("ih1 ~", paste(c("joint", covs), collapse = " + ")))
  one <- function(d) tryCatch({
    m <- geeglm(f, data = d, family = poisson(link = "log"), id = MOMID, corstr = "exchangeable", weights = ipw9)
    list(b = coef(m), V = vcov(m)) }, error = function(e) NULL)
  fits <- mclapply(imp, one, mc.cores = n_cores)
  for (i in which(sapply(fits, is.null))) fits[[i]] <- one(imp[[i]])
  fits <- Filter(Negate(is.null), fits)
  stopifnot(length(fits) == length(imp))
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

# Full model first: must reproduce the current primary estimate
full <- fit_pool("full", red9)
both_full <- with(full[full$exposure == "Both", ], f2(rr, lci, uci))
cat("Full 9-covariate model, both patterns:", both_full, "(target 1.36 (1.04, 1.77))\n")
if (both_full != "1.36 (1.04, 1.77)") {
  if (test_n > 0) message("TEST MODE: not reproduced with ", M, " imputations") else
    stop("The full 9-covariate model does not reproduce aRR 1.36 (1.04, 1.77). Stopping.")
}

crude <- fit_pool("crude", character(0))
add1  <- map(red9, ~ fit_pool(paste0("add_", .x), .x)) %>% set_names(red9)
drop1 <- map(red9, ~ fit_pool(paste0("drop_", .x), setdiff(red9, .x))) %>% set_names(red9)

# ------------------------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------------------------
res <- bind_rows(
  crude %>% mutate(model = "Crude (weighted, no covariates)", type = "crude", covariate = "",
                   base_rr = crude$rr),
  imap_dfr(add1, ~ .x %>% mutate(model = paste("Crude +", lab[[.y]]), type = "add one", covariate = lab[[.y]],
                                 base_rr = crude$rr)),
  imap_dfr(drop1, ~ .x %>% mutate(model = paste("Full -", lab[[.y]]), type = "drop one", covariate = lab[[.y]],
                                  base_rr = full$rr)),
  full %>% mutate(model = "Full (9 covariates)", type = "full", covariate = "", base_rr = crude$rr)) %>%
  mutate(pct_change = 100 * (rr - base_rr) / base_rr,
         change_definition = case_when(type == "crude" ~ "",
                                       type == "add one" ~ "vs weighted crude",
                                       type == "drop one" ~ "vs full",
                                       type == "full" ~ "total change, vs weighted crude"),
         flag_10pct = if_else(type %in% c("add one", "drop one") & abs(pct_change) >= 10, "YES", ""),
         `RR (95% CI)` = f2(rr, lci, uci)) %>%
  transmute(model, type, covariate, exposure, RR = round(rr, 3), lower_ci = round(lci, 3), upper_ci = round(uci, 3),
            `RR (95% CI)`, pct_change = ifelse(type == "crude", NA, round(pct_change, 1)), change_definition,
            flag_10pct, n_cases = 700, n_obs = 40700, n_imps = M)
write_csv(res, file.path(out_dir, "change_in_estimate.csv"))

# Share of the total log-scale change from crude to full, for both patterns
bt <- res %>% filter(exposure == "Both")
tot <- log(full$rr[3]) - log(crude$rr[3])
share <- tibble(covariate = lab[red9],
  add_one_pct = bt$pct_change[bt$type == "add one"], drop_one_pct = bt$pct_change[bt$type == "drop one"],
  add_one_share_of_total = 100 * (log(add1 %>% map_dbl(~ .x$rr[3])) - log(crude$rr[3])) / tot,
  drop_one_share_of_total = 100 * (log(full$rr[3]) - log(drop1 %>% map_dbl(~ .x$rr[3]))) / tot)
print(bt %>% select(model, `RR (95% CI)`, pct_change, flag_10pct), n = Inf, width = Inf)
print(share %>% mutate(across(where(is.numeric), ~ round(.x, 1))), width = Inf)
cat(sprintf("Total change for both: %.1f%% on the RR scale\n", 100 * (full$rr[3] - crude$rr[3]) / crude$rr[3]))
cat("\nDone\n")
