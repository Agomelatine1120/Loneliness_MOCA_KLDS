# Loneliness_MOCA_KLDS

Private reproducible analysis package for the analysis examining whether baseline individual structural covariance similarity moderates the association between longitudinal loneliness change and follow-up cognition.


## Scientific question

Does baseline KL divergence similarity (KLDS) of a prespecified brain edge modify the relationship between change in UCLA Loneliness Scale score and later MoCA performance?

The analysis is edge-wise across 12 prespecified AAL3 edges:

- bilateral hippocampus to right medial OFC, anterior OFC, posterior OFC, and insula;
- bilateral amygdala to the same four targets;
- bilateral dmPFC to the same four targets.

## Model

For each edge:

```text
MOCA_2 ~ MOCA_1 + UCLA_LS_1 + delta_UCLA_LS + KLDS_1
         + delta_UCLA_LS × KLDS_1
         + age_1 + sex + education_1 + FollowupInterval
```

The target coefficient is:

```text
centered(delta_UCLA_LS) × centered(KLDS_1)
```

Because baseline MoCA is included, the outcome is follow-up MoCA conditional on baseline MoCA, corresponding to an ANCOVA/residualized-change interpretation.

## Repository structure

```text
Loneliness_MOCA_KLDS/
├── README.md
├── requirements.txt
├── run_kids.py
├── input/
│   ├── clinical_model3_paired72.tsv
│   ├── baseline_AAL3_12edges_KLDS_paired72.tsv
│   └── AAL3_Model3_12edge_definitions.tsv
├── output/
│   ├── AAL3_KLDS_Model3_12edges_interaction_results.tsv
│   ├── AAL3_KLDS_Model3_12edges_all_coefficients.tsv
│   ├── AAL3_KLDS_Model3_12edges_simple_slopes.tsv
│   ├── AAL3_KLDS_Model3_12edges_paired72_analysis_data.tsv
│   ├── AAL3_KLDS_Model3_12edge_definitions.tsv
│   ├── AAL3_KLDS_Model3_12edges_QC.json
│   └── AAL3_KLDS_Model3_12edges_METHOD.txt
└── docs/
    ├── METHODS_AND_LOGIC.md
    ├── DATA_DICTIONARY.md
    └── RESULTS_SUMMARY.md
```

## Inputs

### `clinical_model3_paired72.tsv`

One row per participant, 72 participants. It contains only the clinical variables required for Model 3 and visit timing.

### `baseline_AAL3_12edges_KLDS_paired72.tsv`

One row per participant, containing `subject_id` plus baseline KLDS for the 12 prespecified edges. LOO, DTW, imaging filenames, and follow-up KLDS are excluded.

### `AAL3_Model3_12edge_definitions.tsv`

The 12 edge keys, seed/target definitions, AAL3 IDs, and AAL3 names.

The included subject identifiers are pseudonymous research codes. Nevertheless, the clinical table is sensitive research data and this repository must remain private unless an authorized data-sharing review approves otherwise.

## Analysis logic

1. Read the clinical, baseline KLDS, and edge-definition tables.
2. Retain exactly three seeds and four right-hemisphere targets.
3. Require exactly 12 unique edges and 72 unique participants.
4. Require the clinical and KLDS participant-ID sets to match exactly.
5. Reshape the data to 72 participants × 12 edges = 864 analysis rows.
6. Check that all Model 3 variables are complete.
7. Within each edge, mean-center baseline MoCA, baseline UCLA-LS, UCLA-LS change, baseline KLDS, age, education, and follow-up interval.
8. Construct the interaction from centered UCLA-LS change and centered baseline KLDS.
9. Fit one ordinary least-squares model per edge.
10. Calculate conventional OLS inference and HC3 robust inference.
11. Calculate unstandardized and standardized coefficients, 95% CIs, partial r, interaction delta-R2, full R2, adjusted R2, and VIF.
12. Check residual normality, heteroskedasticity, and Cook's distance.
13. Estimate simple slopes of UCLA-LS change at baseline KLDS mean -1 SD, mean, and mean +1 SD.
14. Apply BH-FDR across the 12 interaction p-values. A separate 12-test BH-FDR is applied to HC3 p-values as a sensitivity analysis.
15. Write all analysis tables and QC metadata to `output/`.

## Follow-up interval

`FollowupInterval = year_2 - year_1` is included as a continuous covariate in every model.

| Follow-up interval | Participants |
|---:|---:|
| 2 years | 48 |
| 3 years | 18 |
| 4 years | 6 |

## Run the analysis

Create an environment and install the required packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Run from the repository root:

```bash
python3 run_kids.py
```

The script uses repository-relative paths, so no local absolute path editing is required.

## Main outputs

- `AAL3_KLDS_Model3_12edges_interaction_results.tsv`: complete results for all 12 target interactions, including conventional and HC3 BH-FDR.
- `AAL3_KLDS_Model3_12edges_all_coefficients.tsv`: every coefficient in every model, including `FollowupInterval`.
- `AAL3_KLDS_Model3_12edges_simple_slopes.tsv`: 36 simple slopes, three for each edge.
- `AAL3_KLDS_Model3_12edges_paired72_analysis_data.tsv`: final analysis-ready long table.
- `AAL3_KLDS_Model3_12edges_QC.json`: sample, missingness, interval distribution, and result-count verification.

See [`docs/METHODS_AND_LOGIC.md`](docs/METHODS_AND_LOGIC.md) for the complete statistical logic and [`docs/RESULTS_SUMMARY.md`](docs/RESULTS_SUMMARY.md) for the final results.

## Reproducibility notes

- Python 3.9.7 was used for the validated run.
- Sex remains coded as in the source data (0/1); consult the study data dictionary before assigning category labels.
- KLDS is an unsigned similarity measure. Larger values indicate greater similarity between the regional GMV distributions; it is not a signed structural-covariance coefficient.
- The included output files were regenerated from the included input files using `run_kids.py`.
- Conventional OLS results reproduce the original presentation. HC3 inference is essential because several models show heteroskedasticity or influential observations.
