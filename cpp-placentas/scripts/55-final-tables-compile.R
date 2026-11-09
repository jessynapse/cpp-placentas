# ==============================================================================
# 55-final-tables-compile.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Compiles every number for the final manuscript tables (placentas-ppe v3.0_sm)
# from the CURRENT pipeline: 49 (imputation, outcome in the model) -> 50 (IPCW)
# -> 51 (models), 52 (joint reruns), 53 (Boston heterogeneity), 54 (Knol and
# VanderWeele). No models are refitted here. Estimates are read from those
# outputs. Counts, birthweights, E-values and the maternal angioma bounds are
# recomputed from the data and from UNROUNDED estimates.
#
# Inputs: data/weighted_long_m50y.RDS, data/analytic_sample_corrected.RDS,
#         outputs/*_m50y.csv, outputs/rerun*.csv, outputs/interaction_knol_*.csv
# Output: outputs/final_tables_v3.csv
# ==============================================================================

library(tidyverse)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
rd <- function(f) read_csv(file.path(out_dir, f), show_col_types = FALSE)
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)
pct <- function(k, n) sprintf("%.1f%% (%s/%s)", 100 * k / n, fmt_n(k), fmt_n(n))
f2 <- function(r, l, u) sprintf("%.2f (%.2f, %.2f)", r, l, u)
f3 <- function(r, l, u) sprintf("%.3f (%.3f, %.3f)", r, l, u)

# ------------------------------------------------------------------------------
# DATA FOR COUNTS (imputation 1; exposures and outcomes are never imputed)
# ------------------------------------------------------------------------------
w <- readRDS(file.path(in_dir, "weighted_long_m50y.RDS"))
d <- w %>% filter(.imp == 1, !is.na(ih1), exclusion == 0) %>%
  mutate(mvm2n = as.numeric(as.character(mvm2)), ain = as.numeric(as.character(ai)),
         mvm3 = as.character(mvm3), ai3 = as.character(ai3), boston = as.character(site) == "5",
         joint = case_when(mvm2n == 0 & ain == 0 ~ "Neither", mvm2n == 1 & ain == 0 ~ "MVM only",
                           mvm2n == 0 & ain == 1 ~ "AI only", TRUE ~ "Both MVM and AI"))
rm(w); invisible(gc())
stopifnot(nrow(d) == 40700, sum(d$ih1) == 700, sum(d$ih2) == 930)
cnt <- function(x, y = d$ih1) pct(sum(y[x]), sum(x))

# Birthweight by examination status (all 45,268 singleton pregnancies)
a <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS"))
stopifnot(nrow(a) == 45268)
a <- a %>% mutate(grp = case_when(!is.na(ih1) & exclusion == 0 ~ "Examined",
                                  exclusion == 1 ~ "Died before one year",
                                  TRUE ~ "Alive without one-year examination"))
bw <- function(x) sprintf("%.0f (%.0f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
bw_tab <- bind_rows(
  a %>% group_by(grp) %>% summarise(n = n(), n_bw = sum(!is.na(birthweight)), value = bw(birthweight), .groups = "drop"),
  a %>% filter(grp != "Examined") %>% summarise(grp = "All not examined", n = n(),
                                                n_bw = sum(!is.na(birthweight)), value = bw(birthweight)))

# ------------------------------------------------------------------------------
# ESTIMATES FROM THE CURRENT PIPELINE
# ------------------------------------------------------------------------------
ae <- rd("all_estimates_m50y.csv")            # 51
est <- function(analysis, model, term) ae %>% filter(.data$analysis == !!analysis, .data$model == !!model, .data$term == !!term)
e2 <- function(analysis, term) map_chr(c("Crude", "Primary", "Secondary"),
                                        function(m) { x <- est(analysis, m, term); f2(x$rr, x$lci, x$uci) })
syn <- rd("synergy_reri_m50y.csv")            # 51
gs  <- rd("grade_stage_separate_m50y.csv")    # 51
ih2 <- rd("ih2_m50y.csv")                     # 51
bos <- rd("boston_4exp_m50y.csv")             # 51
r1  <- rd("rerun1_ih2_joint.csv")             # 52
r2  <- rd("rerun2_boston_joint.csv")          # 52
kv  <- rd("interaction_knol_strata.csv")      # 54

rows <- list()
add <- function(item, row, ih = "", unadj = "", prim = "", nosite = "", note = "")
  rows[[length(rows) + 1]] <<- tibble(item, row, `IH % (cases/N)` = ih, Unadjusted = unadj,
                                      Primary = prim, `Without study site` = nosite, note)
add3 <- function(item, row, ih, v, note = "") add(item, row, ih, v[1], v[2], v[3], note)

# 1. TABLE 2 ---------------------------------------------------------------------
add("Table 2", "No MVM (reference)", cnt(d$mvm2n == 0), "1.00 (Reference)", "1.00 (Reference)", "1.00 (Reference)")
add3("Table 2", "Any MVM", cnt(d$mvm2n == 1), e2("Any MVM", "mvm2n"))
add("Table 2", "No AI (reference)", cnt(d$ain == 0), "1.00 (Reference)", "1.00 (Reference)", "1.00 (Reference)")
add3("Table 2", "Any AI", cnt(d$ain == 1), e2("Any AI", "ain"))
add("Table 2", "Neither (reference)", cnt(d$joint == "Neither"), "1.00 (Reference)", "1.00 (Reference)", "1.00 (Reference)")
for (j in c("MVM only", "AI only", "Both MVM and AI"))
  add3("Table 2", j, cnt(d$joint == j), e2("Joint", paste0("joint", j)))
for (m in c("RERI", "Attributable proportion", "Ratio of RRs")) {
  x <- syn %>% filter(str_starts(measure, m)); add3("Table 2", m, "", unlist(x[1, 2:4]))
}
for (t in c("mvm2n", "ain")) {
  x <- est(if (t == "mvm2n") "Any MVM" else "Any AI", "Primary", t)
  add("Table 2", paste(if (t == "mvm2n") "Any MVM" else "Any AI", "primary, 3 decimals"), prim = f3(x$rr, x$lci, x$uci))
}

# 2. eTABLE 2 --------------------------------------------------------------------
for (i in seq_len(nrow(kv))) add3("eTable 2", kv$Effect[i], "", unlist(kv[i, 2:4]))

# 3. eTABLE 3A -------------------------------------------------------------------
for (i in seq_len(nrow(gs))) {
  x <- gs[i, ]; ref <- x$level == "None"
  add("eTable 3A", paste(x$analysis, x$level), x$ih_pct,
      if (ref) "1.00 (Reference)" else x$Crude,
      if (ref) "1.00 (Reference)" else x$`Primary: covariates + site + IPW`,
      if (ref) "1.00 (Reference)" else x$`Secondary: covariates + IPW, no site`,
      note = if (x$analysis == "AI stage") "AI stage is the CPP AI_3cat variable (rule undocumented, about 95% match to amnion or umbilical artery = high)" else "")
}
stopifnot(all(table(d$mvm3[d$mvm2n == 1]) > 0))

# 4. eTABLE 3B (suspect IH counted as IH) ----------------------------------------
e2b <- function(lab) { x <- ih2 %>% filter(exposure == lab)
  c(x$rr_ci[x$model == "crude"], x$rr_ci[x$model == "adjusted+site+IPW"], x$rr_ci[x$model == "adjusted+IPW"]) }
add("eTable 3B", "No MVM (reference)", cnt(d$mvm2n == 0, d$ih2), "1.00 (Reference)", "1.00 (Reference)", "1.00 (Reference)")
add3("eTable 3B", "Any MVM", cnt(d$mvm2n == 1, d$ih2), e2b("MVM any vs none"))
add("eTable 3B", "No AI (reference)", cnt(d$ain == 0, d$ih2), "1.00 (Reference)", "1.00 (Reference)", "1.00 (Reference)")
add3("eTable 3B", "Any AI", cnt(d$ain == 1, d$ih2), e2b("AI any vs none"))
for (i in seq_len(nrow(r1))) add3("eTable 3B", r1$`Exposure category`[i], r1$`IH % (cases/N)`[i], unlist(r1[i, 3:5]))

# 5. eTABLE 4 (Boston vs other 11 sites, primary model) --------------------------
for (s in c(TRUE, FALSE)) {
  lab <- if (s) "Boston" else "Other 11 sites"; dd <- d[d$boston == s, ]
  brow <- bos %>% filter(str_starts(stratum, if (s) "Boston" else "Other"))
  add("eTable 4", paste(lab, "stratum"), pct(sum(dd$ih1), nrow(dd)))
  add("eTable 4", paste(lab, "| No MVM (reference)"), cnt(dd$mvm2n == 0, dd$ih1), prim = "1.00 (Reference)")
  add("eTable 4", paste(lab, "| Any MVM"), cnt(dd$mvm2n == 1, dd$ih1), prim = brow$`MVM any vs none`)
  add("eTable 4", paste(lab, "| No AI (reference)"), cnt(dd$ain == 0, dd$ih1), prim = "1.00 (Reference)")
  add("eTable 4", paste(lab, "| Any AI"), cnt(dd$ain == 1, dd$ih1), prim = brow$`AI any vs none`)
  jr <- r2 %>% filter(str_starts(stratum, if (s) "Boston" else "Other"))
  for (i in seq_len(nrow(jr))) add("eTable 4", paste(lab, "|", jr$`Exposure category`[i]), jr$`IH % (cases/N)`[i], prim = jr$`Primary aRR (95% CI)`[i])
}

# 6. eTABLE 5 (E-values from unrounded estimates) --------------------------------
ev <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
ev_ci <- function(l, u) ifelse(l <= 1 & u >= 1, 1, ifelse(l > 1, ev(l), ev(u)))
evx <- list(`Any MVM` = c("Any MVM", "mvm2n"), `Any AI` = c("Any AI", "ain"),
            `MVM only` = c("Joint", "jointMVM only"), `AI only` = c("Joint", "jointAI only"),
            `Both MVM and AI` = c("Joint", "jointBoth MVM and AI"))
evtab <- imap_dfr(evx, function(x, lab) { e <- est(x[1], "Primary", x[2])
  tibble(exposure = lab, rr_ci = f2(e$rr, e$lci, e$uci), ev_point = ev(e$rr), ev_ci = ev_ci(e$lci, e$uci)) })
for (i in seq_len(nrow(evtab))) add("eTable 5", evtab$exposure[i], prim = evtab$rr_ci[i],
  note = sprintf("E-value estimate %.2f, CI limit %.2f", evtab$ev_point[i], evtab$ev_ci[i]))
both <- est("Joint", "Primary", "jointBoth MVM and AI")
bound <- function(rr, ud = 1.64) rr * (ud - 1) / (ud - rr)
add("eTable 5", "Maternal angioma: RR_EU to explain away aRR (RR_UD = 1.64)", note = sprintf("%.2f (aRR %.4f)", bound(both$rr), both$rr))
add("eTable 5", "Maternal angioma: RR_EU to move lower CI limit to 1", note = sprintf("%.2f (lower limit %.4f)", bound(both$lci), both$lci))

# 7. eTABLE 1 (birthweight) ------------------------------------------------------
for (i in seq_len(nrow(bw_tab))) add("eTable 1", paste0("Birthweight, g, mean (SD): ", bw_tab$grp[i], " (n = ", fmt_n(bw_tab$n[i]), ")"),
  note = sprintf("%s, birthweight recorded for %s", bw_tab$value[i], fmt_n(bw_tab$n_bw[i])))

out <- bind_rows(rows)
write_csv(out, file.path(out_dir, "final_tables_v3.csv"))
print(out, n = Inf, width = Inf)
cat("\nE-values (unrounded):\n"); print(evtab)
cat("\nDone\n")
