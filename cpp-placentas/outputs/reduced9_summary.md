# Reduced 9-covariate set (script 71)

## What was done

Three covariates were removed from the primary model: pre-pregnancy BMI, pre-pregnancy diabetes and chronic hypertension. The 9 covariates kept were maternal age, race and ethnicity, education, annual family income, marital status, smoking, parity, infant sex and study site.

The stabilized censoring weights were re-estimated from both placental patterns plus these 9 covariates, with an intercept-only numerator. Fetal deaths, infant deaths and survival without examination stayed as censoring. Everything else was unchanged:

- Same 40,700 examined infants and 700 definite IH cases.
- Same 50 imputations, with no re-imputation.
- Same modified Poisson GEE.
- Same pooling with Rubin's rules.

Before anything new was fit, the 12-covariate pipeline was refit. It reproduced Table 2, including the aRR of 1.37 (1.05, 1.79) for both patterns. All 29 adjusted cells across Table 2 and eTables 2 to 6 matched the current manuscript (v4.6).

The re-estimated weights were well behaved. Their mean was 1.00 and their range was 0.91 to 1.65 in the first imputation, compared with 0.91 to 1.68 for the current weights.

## (a) Does the estimate for both patterns change meaningfully?

No. The aRR for both patterns was 1.36 (1.04, 1.77) with 9 covariates, compared with 1.37 (1.05, 1.79) with 12, a change of 1.0%.

No estimate in any table changed by 10% or more. The largest change in an aRR was 2.5%, for high maternal vascular malperfusion grade, which went from 1.09 (0.76, 1.54) to 1.11 (0.78, 1.58). The largest change of any kind was in the maternal angioma bound in eTable 6, which went from 3.28 to 3.09 (5.8%). This bound is very sensitive to small changes in the aRR, because the aRR is close to the assumed angioma association of 1.64.

## (b) Would any conclusion in Table 2 or the eTables change?

No.

- **Table 2:** both patterns still had the highest risk, with a CI above 1. Each pattern alone was still modestly higher but imprecise: maternal vascular malperfusion only 1.09 (0.90, 1.32) and acute inflammation only 1.05 (0.85, 1.30). The interaction measures stayed imprecise. The RERI was 0.22 (-0.20, 0.63) and the ratio of RRs was 1.18 (0.83, 1.68).
- **eTable 2:** each pattern was still stronger in the presence of the other. For maternal vascular malperfusion, the aRR was 1.29 (0.95, 1.75) with acute inflammation and 1.09 (0.90, 1.32) without it.
- **eTable 3:** there was still no dose-response for either pattern. Low acute inflammation stage still had a CI just above 1, at 1.25 (1.01, 1.53), compared with 1.26 (1.02, 1.55).
- **eTable 4:** the association with both patterns was still stronger at the Boston site, 1.64 (1.16, 2.33), than at the other 11 sites, 1.11 (0.75, 1.65).
- **eTable 5:** counting suspect IH as a case still gave a slightly weaker but more precise estimate for both patterns, 1.28 (1.02, 1.61).
- **eTable 6:** the E-value for both patterns was 2.06 for the estimate and 1.25 for the CI, compared with 2.09 and 1.29.

## (c) Which removal drives the change?

The change comes from removing pre-pregnancy BMI and diabetes, and it is very small. Removing chronic hypertension alone left the aRR for both patterns at exactly 1.37 (1.05, 1.79). Removing BMI and diabetes alone gave 1.36 (1.04, 1.77), the same result as the full 9-covariate set. Every one of these shifts is within rounding error.

## Bottom line

The results do not depend on whether pre-pregnancy BMI, pre-pregnancy diabetes and chronic hypertension are in the model. Keeping the prespecified 12-covariate set as the primary analysis loses nothing. The reduced-set results can be reported as a sensitivity analysis if the committee wants them.

## Files

- `scripts/71-reduced-covariates.R`
- `outputs/reduced9_comparison.csv` (12 vs 9 covariates, and 12 vs each intermediate model)
- `outputs/table2_reduced9.csv`, `outputs/etable2_reduced9.csv` to `outputs/etable6_reduced9.csv`
- `outputs/reduced9_weights.csv`
