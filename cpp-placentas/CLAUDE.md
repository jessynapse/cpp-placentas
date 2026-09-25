# CPP Placental Pathology and Infantile Haemangioma

## Study

Placental pathology and infantile haemangioma (IH) in the Collaborative
Perinatal Project (CPP). Manuscript for the journal *Placenta*. Part of a
PhD dissertation in epidemiology (University of Pennsylvania).

## Sample (corrected pipeline, scripts 07-10, September 2026)

- 45,268 singleton pregnancies: 45,296 in the original sample, minus 26
  multiples (`C10 != 1`) and 2 pregnancies missing AI (complete case on
  both exposures, exposures never imputed)
- 40,700 with an observed IH outcome at one year, all used in models
- 700 definite IH cases (`ih1`)
- 930 when suspect cases count as IH (`ih2`)
- Chronic hypertension 5.8% (2,640), coded `HTN_PREG %in% c(1, 4)`

The original pipeline (01-05, manuscript draft) had 45,296 pregnancies,
40,719 observed, 701 cases, but its adjusted models silently used only
40,683 infants. The corrections changed no estimate by more than 0.01.

## Exposures

- `MVM_SCORE` (continuous, observed range 0 to 8), analysed as `mvm`
- Any MVM (score 1 or more), analysed as `mvm2`
- `AI_DI` (acute inflammation, any vs none), analysed as `ai`
- `AI_3cat` (acute inflammation stage: 0 none, 1 low, 2 high), analysed as `ai3`

## Methods

- Modified Poisson regression using GEE (`geepack::geeglm`, log link)
  with an exchangeable correlation structure, clustered on maternal ID
  (`MOMID`)
- 20 multiply imputed datasets (`mice`, maxit 10). Corrected pipeline
  (08) runs 4 parallel chains of 5 (seeds 12345 to 12348). Gestational
  age and birthweight are not imputation predictors. Models are fit in
  each imputed dataset and pooled with Rubin's rules (manual
  `pool_gee()`, since GEE fits are not `mira` compatible)
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
preeclampsia are also not adjusted for, for the same reason.

**Study site (12 centres) is in the PRIMARY model** (decided with Ellen
Francis, September 2026). Sites differ in population, number of
pathologists and IH ascertainment: any MVM ranges 10% to 80% and definite
IH 0.4% to 3.6% across sites. Site 5 = Boston (41% of cases). Other site
names are unknown (SAS format catalogue not available).

## Model hierarchy and current results (corrected pipeline, `ih1`)

| Model | MVM any | MVM continuous | AI any |
|-------|---------|----------------|--------|
| PRIMARY: 11 covariates + site + IPW | 1.14 (0.96, 1.35) | 1.05 (0.97, 1.15) | 1.14 (0.96, 1.35) |
| SECONDARY: 11 covariates + IPW, no site | 1.25 (1.07, 1.46) | 1.12 (1.03, 1.21) | 1.20 (1.02, 1.42) |
| Crude | 1.22 (1.05, 1.43) | 1.09 (1.01, 1.18) | 1.25 (1.06, 1.47) |

Sensitivity and supplemental: sequential adjustment (13), Boston vs not
Boston (13), suspect-as-present outcome (10), AI stage (12: low 1.32
(1.07, 1.62), high 1.07 (0.84, 1.36), no dose-response).

The original manuscript draft reported the SECONDARY model as primary.
Script 02 descriptive tables and the manuscript counts still need to be
rebuilt from the corrected sample.

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
- `06-site-descriptives.R`: sample, follow-up, IH, exposure and race by
  site, and crude within-site RRs
- `07`-`10`: corrected copies of 01, 03, 04, 05 (changes marked
  `# CHANGED:`). Intermediate files use a `_corrected` suffix.
- `11-compare-original-corrected.R`: original vs corrected, side by side
- `12-ai-stage.R`: AI stage dose-response, IPW rebuilt with `ai3`
- `13-sequential-site.R`: sequential adjustment, primary vs secondary,
  Boston vs not Boston

## Other variables available

- `path2all` has neutrophil infiltration graded 0 to 3 in 7 compartments
  (maternal: chorion PA02_35/37, amnion PA02_34/36. fetal: umbilical
  vein PA02_25, umbilical artery PA02_26, fetal surface vessels PA02_38.
  cord substance PA02_27). `AI_DI` = any of the 7 except cord substance
  (reproduced 100%). `AI_3cat` rule is undocumented (~95% match to
  "high = amnion or umbilical artery involved"), confirm with Alexa.
- Also: bacterial colonies in amnion (PA02_33), decidual neutrophils
  and lymphocytes, macrophages, pathologist knowledge at exam
  (PA101_59, PA02_72), individual MVM lesions, `FVM_DI` (fetal vascular
  malperfusion, not yet used).

## Known issues in the ORIGINAL scripts 01-05 (all fixed in 07-10)

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

5. Gestational age and birthweight were imputation predictors but not
   imputed, and 2 placentas missing AI were treated as complete, so 57
   women kept missing covariates and 36 with observed outcomes were
   silently dropped from adjusted models (N 40,683, not 40,719).
6. Income reference in the code was <=$1,999 (133 women), while the
   manuscript says $4,000-$5,999.

Existing scripts are never edited. Fixes live in new numbered scripts.
