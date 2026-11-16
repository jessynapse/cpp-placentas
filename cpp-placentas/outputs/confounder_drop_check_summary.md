# Confounder drop check (script 70)

## What was done

Six covariates were checked: education, annual family income, marital status, smoking, pre-pregnancy BMI and pre-pregnancy diabetes. Each was dropped from the primary model one at a time. Then the three socioeconomic covariates were dropped together, and then all six. The kept set was maternal age, race and ethnicity, parity, chronic hypertension, infant sex and study site.

Everything else matched Table 2:

- Same 40,700 examined infants and 700 definite IH cases.
- Same four-category exposure, with neither as the reference.
- Same modified Poisson GEE, clustered on mother.
- Same 50 imputations pooled with Rubin's rules.
- Same stabilized censoring weights.

In a last model, the censoring weights were also re-estimated from the reduced covariate set. The full 12-covariate model reproduced Table 2 exactly. A change in the aRR of 10% or more from the full model was flagged.

## Results

| Model | MVM only | AI only | Both | Change for both |
|---|---|---|---|---|
| Crude (weighted, no covariates) | 1.18 (0.98, 1.42) | 1.19 (0.96, 1.47) | 1.56 (1.21, 2.00) | +13.5% |
| Full, 12 covariates | 1.08 (0.89, 1.32) | 1.07 (0.86, 1.32) | 1.37 (1.05, 1.79) | reference |
| Drop education | 1.08 (0.89, 1.32) | 1.06 (0.86, 1.31) | 1.37 (1.05, 1.79) | +0.2% |
| Drop income | 1.08 (0.89, 1.32) | 1.07 (0.86, 1.32) | 1.36 (1.04, 1.78) | -0.8% |
| Drop marital status | 1.08 (0.89, 1.32) | 1.07 (0.86, 1.32) | 1.38 (1.06, 1.79) | +0.3% |
| Drop smoking | 1.08 (0.89, 1.32) | 1.07 (0.87, 1.32) | 1.38 (1.05, 1.79) | +0.2% |
| Drop pre-pregnancy BMI | 1.09 (0.89, 1.32) | 1.06 (0.86, 1.31) | 1.36 (1.05, 1.78) | -0.7% |
| Drop pre-pregnancy diabetes | 1.09 (0.89, 1.32) | 1.06 (0.86, 1.31) | 1.37 (1.05, 1.79) | 0.0% |
| Drop socioeconomic block | 1.08 (0.89, 1.32) | 1.06 (0.86, 1.31) | 1.37 (1.05, 1.78) | -0.5% |
| Reduced set (drop all 6) | 1.09 (0.89, 1.32) | 1.05 (0.85, 1.30) | 1.36 (1.04, 1.77) | -0.9% |
| Reduced set, weights re-estimated | 1.08 (0.89, 1.32) | 1.05 (0.85, 1.30) | 1.35 (1.03, 1.76) | -1.8% |

Values are aRR (95% CI) compared with neither.

## Plain-language answers

**Which covariates move the "both" estimate by 10% or more?** None of the six. Dropping any one of them changed the aRR for both patterns by less than 1%. The largest change for any category was 1.8%. The only flagged changes came from the crude model, which drops all 12 covariates (+13.5% for both and +11.4% for AI only). That gap from crude to adjusted is driven by the covariates that were kept, mainly study site, not by these six.

**Does the reduced set differ materially from the full set?** No. With all six dropped, the aRR for both patterns was 1.36 (1.04, 1.77), compared with 1.37 (1.05, 1.79) in the full model. The single-pattern estimates were also essentially unchanged.

**Does re-estimating the weights from the reduced set change anything?** No. The aRR for both patterns was 1.35 (1.03, 1.76). The re-estimated weights were well behaved, with a mean of 1.00 and a range of 0.91 to 1.51 in the first imputation.

## Why they barely matter

The six covariates are only weakly related to the placental categories. Their largest standardized mean differences compared with neither were all small:

| Covariate | Largest SMD vs neither |
|---|---|
| Education | 0.12 |
| Smoking | 0.12 |
| Income | 0.10 |
| Pre-pregnancy BMI | 0.07 |
| Marital status | 0.05 |
| Pre-pregnancy diabetes | 0.05 |

A covariate can confound only if it is related to both the exposure and the outcome. These six are not related enough to the exposure to shift the estimates, even though some are related to IH.

Crude, unweighted RRs with IH were as follows:

- **Education:** college compared with less than high school, 2.73 (2.09, 3.58).
- **Income:** $15,000 or more compared with $4,000 to $5,999, 3.90 (2.91, 5.23).
- **Marital status:** single compared with married, 0.53 (0.40, 0.69).
- **Smoking:** 20 or more cigarettes a day compared with none, 1.34 (1.11, 1.62).
- **Pre-pregnancy BMI:** per 5 kg/m², 0.76 (0.68, 0.84).
- **Pre-pregnancy diabetes:** 1.86 (1.18, 2.91).

These crude associations may partly reflect differences between study sites in who enrolled and in how completely IH was recorded, which is why site is kept in every model.

Two more descriptive points:

- **Socioeconomic covariates overlap only moderately.** Cramér's V was 0.21 for education and income, 0.15 for education and marital status, and 0.25 for income and marital status. They are not interchangeable, but dropping them together still changed nothing.
- **Pre-pregnancy diabetes is rare.** It was present in 1.5% (600/40,695) of mothers: 1.4% with neither, 2.0% with MVM only, 1.0% with AI only and 1.4% with both. Only 19 infants of mothers with diabetes had IH, 3.2% (19/600), and only 3 of them were in the both category, 7.5% (3/40). With so few cases, diabetes cannot shift the estimates much either way.

## Suggested wording for the committee

The six covariates were included because of prior evidence that they relate to placental pathology, maternal health or socioeconomic position, not because they are established risk factors for IH. Removing them one at a time, as a socioeconomic block, or all together changed the aRR for both patterns by less than 2%, including when the censoring weights were re-estimated. The primary estimate is therefore not sensitive to their inclusion. The main source of confounding in the crude estimate was study site and the other kept covariates.

## Files

- `scripts/70-confounder-drop-check.R`
- `outputs/confounder_drop_check.csv` (one row per model and category)
- `outputs/confounder_descriptives.csv`
