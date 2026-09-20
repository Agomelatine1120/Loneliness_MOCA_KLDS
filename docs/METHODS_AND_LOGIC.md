# Model 3 methods and analysis logic

## Scientific estimand

The interaction asks whether the association between within-person loneliness change and follow-up cognition differs according to the participant's baseline individual SCN similarity.

The model does not establish causal direction. With two clinical time points, it tests prospective moderation while controlling baseline cognition and baseline loneliness.

## Exact model

For each of the 12 edges:

```text
MOCA_2 = beta_0
       + beta_1(MOCA_1 centered)
       + beta_2(UCLA_LS_1 centered)
       + beta_3(delta_UCLA_LS centered)
       + beta_4(KLDS_1 centered)
       + beta_5(delta_UCLA_LS centered × KLDS_1 centered)
       + beta_6(age_1 centered)
       + beta_7(sex_1)
       + beta_8(education_1 centered)
       + beta_9(FollowupInterval centered)
       + error
```

The primary estimand is `beta_5`.

## Why the lower-order terms are required

Both `delta_UCLA_LS` and `KLDS_1` remain in the model whenever their interaction is included. This preserves the hierarchy principle and makes the interaction interpretable as the change in the UCLA-change slope per unit increase in baseline KLDS.

## Why baseline MoCA and baseline UCLA-LS are included

- Adjusting for `MOCA_1` makes `MOCA_2` a residualized-change/ANCOVA outcome.
- Adjusting for `UCLA_LS_1` separates loneliness change from the participant's starting loneliness level.
- This is generally more stable than modeling a raw MoCA change score alone.

## Why FollowupInterval is included

The observed follow-up interval is not constant: 48 participants were followed for 2 years, 18 for 3 years, and 6 for 4 years. Centered `FollowupInterval` is therefore included in every edge model.

## Edge-wise family

The hypothesis family contains exactly 12 prespecified interactions:

```text
3 bilateral seeds × 4 right targets = 12 tests
```

BH-FDR is applied across these 12 conventional OLS interaction p-values. HC3 p-values receive a separate 12-test BH-FDR sensitivity correction.

## Effect-size and model outputs

For the interaction, the script reports:

- raw coefficient and standardized coefficient;
- standard error, 95% CI, t statistic, and p-value;
- partial correlation calculated from t and residual degrees of freedom;
- delta-R2 compared with the same model without the interaction;
- full-model R2 and adjusted R2.

## Simple slopes

For every edge, the association between UCLA-LS change and follow-up MoCA is estimated at:

- baseline KLDS mean -1 SD;
- baseline KLDS mean;
- baseline KLDS mean +1 SD.

A positive interaction means that the UCLA-change slope becomes more positive as baseline KLDS increases. The sign alone does not prove vulnerability or resilience; the simple slopes determine whether the pattern is buffering, amplifying, or cross-over.

## Diagnostics and sensitivity analysis

The script calculates:

- predictor VIF;
- Shapiro-Wilk residual p-value;
- Breusch-Pagan heteroskedasticity p-value;
- maximum Cook's distance;
- number of observations above `4/n` Cook's-distance threshold.

Because heteroskedasticity and influential observations are present in several models, HC3 inference is treated as an essential sensitivity analysis rather than an optional supplement.

## Interpretation limits

- KLDS is an unsigned similarity measure and is not a signed covariance coefficient.
- An interaction is not evidence that loneliness change causes cognitive change.
- A positive high-KLDS simple slope should not be stated as evidence that increasing loneliness improves cognition.
- MoCA ceiling effects, regression to the mean, and influential cases are plausible competing explanations.
- Findings that survive both conventional and HC3 FDR provide the strongest evidence.

