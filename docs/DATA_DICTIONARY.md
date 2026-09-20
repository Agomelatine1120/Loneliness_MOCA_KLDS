# Data dictionary

## Clinical input

File: `input/clinical_model3_paired72.tsv`

| Variable | Definition | Analysis use |
|---|---|---|
| `subject_id` | Pseudonymous participant identifier | Clinical-imaging join key |
| `year_1` | Baseline assessment year | Derivation/QC |
| `year_2` | Follow-up assessment year | Derivation/QC |
| `FollowupInterval` | `year_2 - year_1`, in years | Covariate |
| `age_1` | Age at baseline | Covariate |
| `sex_1` | Sex in source 0/1 coding | Covariate |
| `education_1` | Baseline years of education | Covariate |
| `UCLA_LS_1` | Baseline UCLA Loneliness Scale | Baseline adjustment |
| `UCLA_LS_2` | Follow-up UCLA Loneliness Scale | Derivation/QC |
| `delta_UCLA_LS` | `UCLA_LS_2 - UCLA_LS_1` | Longitudinal exposure |
| `MOCA_1` | Baseline Montreal Cognitive Assessment | Baseline outcome adjustment |
| `MOCA_2` | Follow-up Montreal Cognitive Assessment | Dependent variable |

## Baseline KLDS input

File: `input/baseline_AAL3_12edges_KLDS_paired72.tsv`

The table contains `subject_id` and 12 columns named:

```text
<edge_key>__KLDS
```

KLDS is calculated upstream as `exp[-(KL(P||Q) + KL(Q||P))]`. Higher KLDS indicates greater similarity between the two regional gray-matter distributions. The bilateral-seed value is the arithmetic mean of the left-seed-to-target and right-seed-to-target KLDS values.

## Edge-definition input

File: `input/AAL3_Model3_12edge_definitions.tsv`

| Field | Meaning |
|---|---|
| `edge_key` | Machine-readable edge identifier |
| `edge` | Human-readable edge description |
| `seed` | Seed family: hippocampus, amygdala, or dmPFC |
| `target` | Right target region |
| `seed_left_AAL3_ID` | Left seed AAL3 label ID |
| `seed_right_AAL3_ID` | Right seed AAL3 label ID |
| `target_AAL3_ID` | Target AAL3 label ID |
| `seed_left_AAL3_name` | Left seed AAL3 label name |
| `seed_right_AAL3_name` | Right seed AAL3 label name |
| `target_AAL3_name` | Target AAL3 label name |

## Derived regression variables

All continuous predictors are centered separately within each edge-specific data set:

```text
moca1_c          = MOCA_1 - mean(MOCA_1)
ucla1_c          = UCLA_LS_1 - mean(UCLA_LS_1)
delta_ucla_c     = delta_UCLA_LS - mean(delta_UCLA_LS)
baseline_klds_c  = KLDS_1 - mean(KLDS_1)
age1_c           = age_1 - mean(age_1)
education1_c     = education_1 - mean(education_1)
interval_c       = FollowupInterval - mean(FollowupInterval)
interaction      = delta_ucla_c × baseline_klds_c
```

Sex is not centered and retains the source 0/1 coding.

