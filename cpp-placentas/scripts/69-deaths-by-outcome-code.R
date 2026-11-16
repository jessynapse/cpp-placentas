# ==============================================================================
# 69-deaths-by-outcome-code.R
# CPP Placental Pathology -> Infantile Haemangioma
# November 2026
#
# Are stillbirths among the 1,834 infants counted as "died before 12 months"?
# 07 defines exclusion = 1 when C1092 ("outcome") is not 0 and not 51-58.
# The C1092 value labels are in a SAS format catalog that is not available,
# so codes are identified from other variables: 1-minute Apgar (C555, C560),
# first cry after delivery (C1122), infant death (C1091), outcome code (C1094),
# gestational age and birthweight. A fetal death has no Apgar and no cry.
#
# Output: outputs/deaths_by_outcome_code.csv (aggregate counts only)
# ==============================================================================

library(tidyverse)
library(haven)

in_dir  <- if (dir.exists("/Users/wongjj/Downloads")) "/Users/wongjj/Downloads" else "data"
out_dir <- if (basename(getwd()) == "scripts") "../outputs" else "outputs"
cn <- function(x) { x <- as.character(x); x[x %in% c(".u", ".U", "u", "U")] <- NA; suppressWarnings(as.numeric(x)) }

a <- readRDS(file.path(in_dir, "analytic_sample_corrected.RDS")) %>%
  select(MOMID, PREGID, CHILDID, ih1, exclusion, birthweight, gest_age)
stopifnot(nrow(a) == 45268, sum(a$exclusion == 1) == 1834)
n <- read_sas(file.path(in_dir, "ncppbasa.sas7bdat"),
              col_select = c(MOMID, PREGID, CHILDID, C1091, C1092, C1094, C555, C560, C1122))
x <- a %>% left_join(n, by = c("MOMID", "PREGID", "CHILDID")) %>%
  mutate(across(c(C1091, C1092, C1094, C555, C560, C1122), cn))

out <- x %>% filter(exclusion == 1) %>%
  group_by(outcome_C1092 = C1092) %>%
  summarise(n = n(), outcome_code_C1094 = paste(sort(unique(C1094)), collapse = ","),
            apgar1_recorded = sum(!is.na(C560)), first_cry_yes = sum(C1122 == 1, na.rm = TRUE),
            median_ga_wk = median(gest_age, na.rm = TRUE), median_bw_g = median(birthweight, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(classification = if_else(apgar1_recorded == 0 & first_cry_yes == 0,
                                   "Fetal death (no Apgar, no cry)", "Live birth, died before 12 months"))
print(out, n = Inf, width = Inf)
cat("\nTotals:\n"); print(out %>% group_by(classification) %>% summarise(n = sum(n)))
dir.create(out_dir, showWarnings = FALSE)
write_csv(out, file.path(out_dir, "deaths_by_outcome_code.csv"))
