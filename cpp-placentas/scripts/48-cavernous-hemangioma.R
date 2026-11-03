# ==============================================================================
# 48-cavernous-hemangioma.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# How often the 12-month exam (PED-12) recorded "cavernous hemangioma"
# (C984) alongside "strawberry hemangioma" (C983, the study outcome), and
# how the newborn nursery record (PED-8: C737 strawberry/portwine, C738
# cavernous) relates to the 12-month diagnosis. Analytic sample (40,700).
#
# Coding of the PED-12 items follows 07: 1 = definite, 2 = suspect,
# 0 = absent (values printed below to confirm).
#
# Inputs: data/analytic_sample_corrected.RDS (07), data/ncppbasa.sas7bdat
# Output: outputs/cavernous_hemangioma.csv
# ==============================================================================

library(tidyverse)
library(haven)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

d <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  filter(!is.na(ih1), exclusion == 0) %>% select(MOMID, PREGID, CHILDID, ih1, mvm2, ai)
stopifnot(nrow(d) == 40700, sum(d$ih1) == 700)
ncpp <- read_sas(file.path(in_dir, "ncppbasa.sas7bdat"),
                 col_select = c(MOMID, PREGID, CHILDID, C982, C983, C984, C737, C738))
d <- d %>% left_join(ncpp, by = c("MOMID", "PREGID", "CHILDID")) %>%
  mutate(across(c(C982, C983, C984, C737, C738), num))

cat("Raw values, analytic sample:\n")
for (v in c("C983", "C984", "C982", "C737", "C738")) { cat(v, ": "); print(table(d[[v]], useNA = "ifany")) }

lab3 <- function(x) case_when(x == 1 ~ "Definite", x == 2 ~ "Suspect", x == 0 ~ "Absent", TRUE ~ "Missing")
cat("\nStrawberry (C983, rows) by cavernous (C984, columns):\n")
print(table(Strawberry = lab3(d$C983), Cavernous = lab3(d$C984)))

cav_def <- d$C984 == 1 & !is.na(d$C984)
cav_any <- d$C984 %in% c(1, 2)
out <- tibble(
  measure = c("Definite strawberry IH (study outcome)",
              "Definite cavernous haemangioma",
              "Definite cavernous, also definite strawberry",
              "Definite cavernous, not definite strawberry",
              "Definite or suspect cavernous",
              "Port-wine stain, definite (not IH)"),
  n = c(sum(d$ih1), sum(cav_def), sum(cav_def & d$ih1 == 1), sum(cav_def & d$ih1 == 0),
        sum(cav_any), sum(d$C982 == 1, na.rm = TRUE))) %>%
  mutate(pct_of_40700 = round(100 * n / nrow(d), 2))
print(out, width = Inf)

cat("\nNewborn record (PED-8 strawberry/portwine, C737) by 12-month definite IH:\n")
print(table(Newborn = lab3(d$C737), IH12 = d$ih1, useNA = "ifany"))

dir.create(out_dir, showWarnings = FALSE)
write_csv(out, file.path(out_dir, "cavernous_hemangioma.csv"))
cat("Saved:", file.path(out_dir, "cavernous_hemangioma.csv"), "\n")
