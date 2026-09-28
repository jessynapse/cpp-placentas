# ==============================================================================
# 27-site-table-percentages.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Site table for Stefanie Hinkle: IH cases per site with percentages and
# the site terms from the primary model (from 22).
# Site names: 5 = Boston is confirmed (05-outcomemodels.R). Other names are
# from Jessica's mapping and are consistent with each site's demographics,
# but should be confirmed against the CPP codebook or with Alexa Freedman.
#
# Inputs: analytic_sample_corrected.RDS, outputs/site_coefficients.csv
# Output: outputs/site_table_percentages.csv
# ==============================================================================

library(tidyverse)

in_dir  <- "/Users/wongjj/Downloads"
out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")

site_levels <- c("5", "10", "15", "31", "37", "45", "50", "55", "60", "66", "71", "82")
site_names  <- c("5" = "Boston", "10" = "Buffalo", "15" = "New Orleans",
                 "31" = "New York (Columbia)", "37" = "Baltimore",
                 "45" = "Virginia", "50" = "Minnesota",
                 "55" = "New York Medical College", "60" = "Oregon",
                 "66" = "Pennsylvania", "71" = "Providence", "82" = "Tennessee")
fmt <- function(x, n) sprintf("%.1f%% (%s/%s)", 100 * x / n,
                              format(x, big.mark = ",", trim = TRUE),
                              format(n, big.mark = ",", trim = TRUE))

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  filter(!is.na(ih1), exclusion == 0) %>%
  mutate(site = factor(as.character(site), levels = site_levels))

tab <- d %>% group_by(site) %>%
  summarise(
    `IH cases` = fmt(sum(ih1 == 1), n()),
    .groups = "drop") %>%
  mutate(site = as.character(site))

coefs <- read_csv(file.path(out_dir, "site_coefficients.csv"), show_col_types = FALSE) %>%
  mutate(site = sub(" .*", "", site)) %>%
  select(site, `Site RR (MVM model)` = `Primary model, MVM any`,
         `Site RR (AI model)` = `Primary model, AI any`)

out <- tab %>% left_join(coefs, by = "site") %>%
  mutate(site = paste(site, site_names[site])) %>%
  rename(Site = site) %>%
  mutate(across(starts_with("Site RR"), ~ ifelse(.x == "1.00 (ref)", "1 REF", .x)))

print(out, width = Inf)
write_csv(out, file.path(out_dir, "site_table_percentages.csv"))
cat("\nSaved: outputs/site_table_percentages.csv\n")
