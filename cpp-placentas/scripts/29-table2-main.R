# ==============================================================================
# 29-table2-main.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Main manuscript Table 2 (definite IH, 40,700 infants, 700 cases), per
# Ellen Francis (September 2026): single-exposure models for MVM and AI
# plus the joint MVM x AI exposure moved from the supplement into the main
# table. Columns: crude, primary (11 covariates + site + IPW), secondary
# (11 covariates + IPW, no site).
#
# Assembled from existing results (no new models):
#   outputs/four_exposures_results.csv   (14)
#   outputs/corrected_group_counts.csv   (10)
#   outputs/joint_models.csv, joint_counts.csv (16)
#
# Output: outputs/table2_main.csv. No p-values.
# ==============================================================================

library(tidyverse)

out_dir <- file.path(normalizePath(file.path(getwd(), "..")), "outputs")
if (basename(getwd()) != "scripts") out_dir <- file.path(getwd(), "outputs")
rd <- function(f) read_csv(file.path(out_dir, f), show_col_types = FALSE)

# single exposures
single <- rd("four_exposures_results.csv") %>%
  select(exposure, model, rr_ci) %>%
  pivot_wider(names_from = model, values_from = rr_ci) %>%
  transmute(label = recode(exposure,
                           "MVM any vs none" = "Any MVM vs none",
                           "MVM score, per 1 point (0 to 8)" = "MVM score, per 1 point",
                           "AI any vs none" = "Any AI vs none",
                           "AI score, per 1 compartment (0 to 7)" = "AI score, per 1 compartment"),
            Crude, Primary = `Primary: + site + IPW`, Secondary = `Secondary: no site + IPW`)

gc <- rd("corrected_group_counts.csv") %>% filter(outcome == "ih1", model == "crude")
ih_of <- function(v, lvl) gc$ih_pct[gc$variable == v & gc$level == lvl]
single <- single %>%
  mutate(IH = case_when(
    label == "Any MVM vs none" ~ paste0(ih_of("mvm2_num", "1"), " vs ", ih_of("mvm2_num", "0")),
    label == "Any AI vs none"  ~ paste0(ih_of("ai_num", "1"), " vs ", ih_of("ai_num", "0")),
    TRUE ~ "1.7% (700/40,700)"))

# joint exposure
jc <- rd("joint_counts.csv") %>% filter(table == "IH by joint group")
joint <- rd("joint_models.csv") %>%
  filter(analysis == "Joint 4-level exposure") %>%
  mutate(model = case_when(grepl("^Crude", model) ~ "Crude",
                           grepl("^Primary", model) ~ "Primary",
                           TRUE ~ "Secondary"),
         label = recode(contrast, "MVM only vs neither" = "MVM only",
                        "AI only vs neither" = "AI only", "Both vs neither" = "Both MVM and AI")) %>%
  select(label, model, rr_ci) %>%
  pivot_wider(names_from = model, values_from = rr_ci)
joint <- bind_rows(tibble(label = "Neither (reference)", Crude = "1 REF",
                          Primary = "1 REF", Secondary = "1 REF"), joint) %>%
  mutate(IH = jc$value[match(sub(" \\(reference\\)", "", label),
                             sub("Both MVM and AI", "Both MVM and AI", jc$measure))])

table2 <- bind_rows(
  tibble(label = "Single exposures"), single,
  tibble(label = "Joint exposure"), joint) %>%
  select(Exposure = label, `IH, % (cases/N)` = IH, Crude,
         `Primary (covariates + site)` = Primary,
         `Secondary (covariates, no site)` = Secondary) %>%
  mutate(across(everything(), ~ replace_na(.x, "")))

print(table2, n = Inf, width = Inf)
write_csv(table2, file.path(out_dir, "table2_main.csv"))
cat("\nSaved: outputs/table2_main.csv\n")
