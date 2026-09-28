# ==============================================================================
# 27-site-table-percentages.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Site table for Stefanie Hinkle: IH cases per site with percentages, IH
# by MVM and by AI within each site (to show sparse cells), and the site
# terms from the primary model (from 22).
# Percentages are IH risk within each group: cases / infants in group.
#
# Inputs: analytic_sample_corrected.RDS, outputs/site_coefficients.csv
# Output: outputs/site_table_percentages.csv
# ==============================================================================

library(tidyverse)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")

site_levels <- c("5", "10", "15", "31", "37", "45", "50", "55", "60", "66", "71", "82")
fmt <- function(x, n) sprintf("%.1f%% (%s/%s)", 100 * x / n,
                              format(x, big.mark = ",", trim = TRUE),
                              format(n, big.mark = ",", trim = TRUE))

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  filter(!is.na(ih1), exclusion == 0) %>%
  mutate(site = factor(as.character(site), levels = site_levels))

tab <- d %>% group_by(site) %>%
  summarise(
    N = format(n(), big.mark = ","),
    `IH, all` = fmt(sum(ih1 == 1), n()),
    `IH, MVM` = fmt(sum(ih1 == 1 & mvm2 == 1), sum(mvm2 == 1)),
    `IH, no MVM` = fmt(sum(ih1 == 1 & mvm2 == 0), sum(mvm2 == 0)),
    `IH, AI` = fmt(sum(ih1 == 1 & ai == 1), sum(ai == 1)),
    `IH, no AI` = fmt(sum(ih1 == 1 & ai == 0), sum(ai == 0)),
    .groups = "drop") %>%
  mutate(site = as.character(site))

coefs <- read_csv(file.path(out_dir, "site_coefficients.csv"), show_col_types = FALSE) %>%
  mutate(site = sub(" .*", "", site)) %>%
  select(site, `Site RR (MVM model)` = `Primary model, MVM any`,
         `Site RR (AI model)` = `Primary model, AI any`)

out <- tab %>% left_join(coefs, by = "site") %>%
  mutate(site = ifelse(site == "5", "5 Boston", site))

print(out, width = Inf)
write_csv(out, file.path(out_dir, "site_table_percentages.csv"))
cat("\nSaved: outputs/site_table_percentages.csv\n")
