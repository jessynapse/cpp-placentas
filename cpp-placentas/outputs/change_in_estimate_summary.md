# Which covariates drive the change in RR? (script 73)

## What was done

This analysis used the 9-covariate primary model:

- 40,700 examined infants with 700 definite IH cases.
- 50 imputations pooled with Rubin's rules.
- Modified Poisson GEE clustered on mother.
- The four-category placental injury pattern as the exposure, with neither as the reference.

Every model used the same stabilized censoring weights as the primary model, estimated from both placental patterns plus the 9 covariates. Any change in the estimates is therefore due to the covariates, not to the weighting.

Four kinds of model were fit:

- **Crude:** no covariates.
- **Add one:** crude plus each covariate alone.
- **Drop one:** the full model minus each covariate.
- **Full:** all 9 covariates.

The full model reproduced the current primary estimate for both patterns, aRR 1.36 (1.04, 1.77).

## Results for both patterns vs neither

| Covariate | Crude + covariate, RR (95% CI) | Add-one change | Full minus covariate, RR (95% CI) | Drop-one change |
|---|---|---|---|---|
| Maternal age | 1.54 (1.20, 1.98) | -0.9% | 1.36 (1.04, 1.78) | +0.2% |
| Race and ethnicity | 1.49 (1.16, 1.92) | -3.9% | 1.37 (1.05, 1.79) | +0.9% |
| Education | 1.61 (1.25, 2.07) | +3.8% | 1.36 (1.04, 1.77) | +0.1% |
| Annual family income | 1.58 (1.23, 2.03) | +1.8% | 1.35 (1.03, 1.76) | -0.7% |
| Marital status | 1.56 (1.21, 2.00) | +0.4% | 1.36 (1.05, 1.78) | +0.3% |
| Smoking | 1.53 (1.19, 1.97) | -1.3% | 1.36 (1.04, 1.78) | +0.3% |
| Parity | 1.59 (1.23, 2.04) | +2.4% | 1.36 (1.04, 1.77) | +0.1% |
| Infant sex | 1.58 (1.23, 2.04) | +2.1% | 1.33 (1.02, 1.73) | -2.3% |
| **Study site** | **1.32 (1.01, 1.73)** | **-14.7%** | **1.56 (1.22, 2.01)** | **+15.0%** |

- The weighted crude RR for both patterns was 1.55 (1.21, 2.00).
- The full 9-covariate aRR was 1.36 (1.04, 1.77).
- Add-one changes are relative to the weighted crude RR. Drop-one changes are relative to the full aRR.

## Plain-language answers

**(a) Total change from crude to adjusted.** For both patterns, adjustment lowered the RR from 1.55 (1.21, 2.00) to 1.36 (1.04, 1.77), a change of -12.4%. The single-pattern estimates also fell:

- Maternal vascular malperfusion only: from 1.18 (0.98, 1.42) to 1.09 (0.90, 1.32), -7.8%.
- Acute inflammation only: from 1.19 (0.96, 1.47) to 1.05 (0.85, 1.30), -11.3%.

**(b) Which one or two covariates account for most of the change.** Study site accounts for essentially all of it, and both approaches agree.

- **Add one:** adding study site alone to the crude model lowered the RR for both patterns by 14.7%. That is more than the full adjustment did.
- **Drop one:** removing study site from the full model raised the RR by 15.0%, back to the crude value.

The CPP sites differed widely in how often placentas showed these patterns and in how often IH was recorded, so site acts as a strong confounder.

No other covariate comes close. Race and ethnicity is the next largest on the add-one side, at -3.9%. Most of that is already captured by site, because dropping race from the full model changed the estimate by only 0.9%. Infant sex is the next largest on the drop-one side, at -2.3%, which is small. Several covariates pushed the crude estimate slightly upward when added alone, by 2% to 4%: education, income, parity and infant sex. Once site was in the model, removing any of them changed almost nothing.

**(c) Does any single covariate move the estimate for both patterns by 10% or more?** Only study site, at -14.7% when added and +15.0% when dropped. Every other covariate moved the estimate by less than 4% in either direction.

## Answer for the committee chair

Yes. One covariate, study site, is responsible for nearly all of the change in RR. Adjusting for site alone moves the RR for both patterns from 1.55 to 1.32, and removing site from the full model moves it back to 1.56. No other covariate changes the estimate by more than 4% in either analysis. This supports keeping study site in the primary model, as decided with Ellen Francis.

## Files

- `scripts/73-change-in-estimate.R`
- `outputs/change_in_estimate.csv` (one row per model and exposure category)
