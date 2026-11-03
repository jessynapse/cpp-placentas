# ==============================================================================
# 47-age-at-ih-exam.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Age at the 12-month paediatric examination (PED-12), at which IH (C983,
# "strawberry hemangioma (ped-12)") was recorded. C809 = "age (wks) at
# 1-year neurological exam", recorded at the same PED-12 visit. Summarised
# for the 40,700 infants in the analytic sample, overall and by definite IH.
#
# Inputs: data/analytic_sample_corrected.RDS (07), data/ncppbasa.sas7bdat
# Output: outputs/age_at_ih_exam.csv
# ==============================================================================

library(tidyverse)
library(haven)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  filter(!is.na(ih1), exclusion == 0) %>%
  select(MOMID, PREGID, CHILDID, ih1)
stopifnot(nrow(d) == 40700, sum(d$ih1) == 700)

ncpp <- read_sas(file.path(in_dir, "ncppbasa.sas7bdat"),
                 col_select = c(MOMID, PREGID, CHILDID, C808, C809))
d <- d %>% left_join(ncpp, by = c("MOMID", "PREGID", "CHILDID")) %>%
  mutate(age_wk = suppressWarnings(as.numeric(as.character(C809))))

summ <- function(x, grp) {
  x <- x[!is.na(x)]
  tibble(group = grp, n_with_age = length(x),
         mean_wk = round(mean(x), 1), sd_wk = round(sd(x), 1),
         median_wk = median(x), q1_wk = quantile(x, .25), q3_wk = quantile(x, .75),
         min_wk = min(x), max_wk = max(x),
         pct_48_to_56wk = round(100 * mean(x >= 48 & x <= 56), 1),
         pct_under_48wk = round(100 * mean(x < 48), 1),
         pct_over_56wk  = round(100 * mean(x > 56), 1),
         mean_months = round(mean(x) * 7 / 30.4375, 1))
}
out <- bind_rows(summ(d$age_wk, "All examined (40,700)"),
                 summ(d$age_wk[d$ih1 == 1], "Definite IH"),
                 summ(d$age_wk[d$ih1 == 0], "No IH")) %>%
  mutate(n_missing_age = c(sum(is.na(d$age_wk)), sum(is.na(d$age_wk[d$ih1 == 1])),
                           sum(is.na(d$age_wk[d$ih1 == 0]))))
print(out, width = Inf)
cat("\nDistribution of age (weeks), all examined:\n")
print(quantile(d$age_wk, c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1), na.rm = TRUE))

dir.create(out_dir, showWarnings = FALSE)
write_csv(out, file.path(out_dir, "age_at_ih_exam.csv"))
cat("Saved:", file.path(out_dir, "age_at_ih_exam.csv"), "\n")
