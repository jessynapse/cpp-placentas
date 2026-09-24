# CPP Placental Pathology and Infantile Haemangioma

## Study

Placental pathology and infantile haemangioma (IH) in the Collaborative
Perinatal Project (CPP). Manuscript for the journal *Placenta*. Part of a
PhD dissertation in epidemiology (University of Pennsylvania).

## Sample

- 45,296 singleton pregnancies (but see Known issues: 26 are twins)
- 40,719 with an observed IH outcome at one year
- 701 definite IH cases (`ih1`)
- 931 when suspect cases count as IH (`ih2`)

## Exposures

- `MVM_SCORE` (continuous, observed range 0 to 8), analysed as `mvm`
- Any MVM (score 1 or more), analysed as `mvm2`
- `AI_DI` (acute inflammation, any vs none), analysed as `ai`
- `AI_3cat` (acute inflammation stage: 0 none, 1 low, 2 high), analysed as `ai3`

## Methods

- Modified Poisson regression using GEE (`geepack::geeglm`, log link)
  with an exchangeable correlation structure, clustered on maternal ID
  (`MOMID`)
- 20 multiply imputed datasets (`mice`, maxit 10, seed 12345), pooled
  with Rubin's rules (manual `pool_gee()` in script 05, since GEE fits
  are not `mira` compatible)
- Stabilised inverse probability of censoring weights. Numerator is
  intercept only. Denominator includes `mvm2`, `ai`, all 11 covariates,
  and `site`. In the code, "censored" (`ih_obs == 0`) means died before
  one year OR survived without an observed IH outcome, so the weights
  cover loss to follow-up as well as non-survival, even though the
  script header says LTFU is not weighted.
- Weights are computed within each imputed dataset. Outcome models are
  fit on `ih_obs == 1` only.

## Covariates (11)

Maternal age, pre-pregnancy BMI, race and ethnicity, education, income,
marital status, smoking, parity, pre-pregnancy diabetes, chronic
hypertension, infant sex.

Gestational age and birthweight are NOT adjusted for in primary models
because they are potential mediators. Gestational hypertension and
preeclampsia are also not adjusted for, for the same reason. Site is
not in primary models (it captures pathologist variability) but is a
sensitivity analysis in script 05.

## Primary results to reproduce (adjusted + IPW, `ih1`)

- MVM any: aRR 1.25 (1.07, 1.46)
- MVM continuous: aRR 1.12 (1.03, 1.21)
- AI any: aRR 1.20 (1.02, 1.42)

## Rules for all output

- Never report p-values
- Format estimates as aRR 1.25 (1.07, 1.46), two decimals
- Report percentages before counts, e.g. 2.0% (225/11,334)
- Never overwrite or edit existing scripts or raw data. Put new work in
  new numbered scripts (06, 07, and so on)
- Save every results table as a CSV in an `outputs/` folder
- Always report the number of cases and the number of observations for
  every model
- Never commit data files. `data/` and all `.csv`, `.rds` and
  `.sas7bdat` files are git-ignored for a reason.
- In any prose written for the user, never use em-dashes or semicolons

## Code environment

- Language: R
- Packages: `haven`, `tidyverse`, `labelled` (01), `table1`, `tableone`
  (02), `mice`, `VIM` (03), `mice` (04), `geepack` (05)
- Every existing script starts with `setwd("/Users/wongjj/Downloads")`
  and reads and writes files there. On the user's Mac, that is where the
  data live.
- In Claude Code cloud sessions, raw data are in `data/` (git-ignored)
  and must be re-uploaded each session. Large files (for example the
  233 MB `ncppbasa.sas7bdat`) are uploaded zipped.

## Data files

Raw inputs (all read by script 01):

| File | Contents |
|------|----------|
| `path1all.sas7bdat` | Placental pathology, part 1 |
| `path2all.sas7bdat` | Placental pathology, part 2 |
| `ncppbasa.sas7bdat` | Main CPP file (outcome, covariates, plurality `C10`, site) |
| `cpp_diabetes_20231026.sas7bdat` | `DM` by `MOMID` and `PREGID` |
| `cpp_htn_20231031.sas7bdat` | `HTN_PREG`: 0 normotensive, 1 chronic, 2 gestational HTN, 3 preeclampsia, 4 superimposed preeclampsia |
| `CPP_dataset_021226.csv` | MVM and AI scores, keyed by `NINDB_NUMBER` (= `PREGID`) |

Intermediate files created by the pipeline: `analytic_sample.RDS` (01),
`imputed_long.csv` and `.RDS` (03), `weighted_long.RDS` (04),
`outcome_results.csv` (05).

## What each script does

- `01-load-clean.R`: loads and merges raw files, defines IH outcomes
  (`ih1` suspect = absent, `ih2` suspect = present, both NA if died or
  excluded by `C1092`), excludes twins and missing MVM, derives
  covariates (parity via `na_tag()` so `.P` = 0 and `.U` = NA), sex
  specific SGA, and saves `analytic_sample.RDS`
- `02-tableone.R`: Table 1 by MVM and by AI, supplemental table by IH
  status (`supp_table_by_ih.csv`), cross-tabs, and a differential
  exclusion check
- `03-imputation.R`: MICE on covariates only (outcomes never imputed),
  saves the long imputed dataset
- `04-ipw.R`: stabilised weights for having an observed IH outcome,
  within each imputation, plus weight diagnostics
- `05-outcomemodels.R`: crude, primary (adjusted and adjusted + IPW),
  site sensitivity, `ih2` sensitivity, and supplemental models (3-level
  MVM and AI, MVM x AI, chorangioma), pooled with Rubin's rules

## Known issues (found September 2026, not yet fixed)

1. **26 twins remain in the analytic sample.** Script 01 identifies
   twins only as `PREGID` appearing more than once after merging. When a
   co-twin has no pathology record, the other twin passes as a
   singleton. One of the 26 is a definite IH case. Fix: also require
   `C10 == 1`.
2. **Chronic hypertension is under-counted.** Script 01 codes
   `chronic_htn = HTN_PREG == 1`, but level 4 (superimposed
   preeclampsia) also means chronic hypertension. This puts 888 women
   with chronic HTN in the reference group. Fix:
   `HTN_PREG %in% c(1, 4)`.
3. The header note in script 01 says "31 singleton conflicts" in the MVM
   file. In the February 2026 file, all 194 pregnancy IDs with
   conflicting MVM records are twins or triplets, and none are
   singletons.
4. The chronic HTN variable label still says `HTN_PREG==1`.

Because existing scripts must not be edited, fixes for these go in new
numbered scripts.
