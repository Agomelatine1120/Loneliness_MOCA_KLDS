#!/usr/bin/env python3
"""Final Model 3-only analysis for 12 focused AAL3 KLDS edges.

Model:
    MOCA_2 ~ MOCA_1 + UCLA_LS_1 + delta_UCLA_LS + KLDS_1
             + delta_UCLA_LS:KLDS_1
             + age_1 + sex + education_1 + FollowupInterval

The primary target is the centered delta_UCLA_LS x centered KLDS_1
interaction. BH-FDR is applied across the 12 prespecified edges.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor


PROJECT_DIR = Path(__file__).resolve().parent
INPUT_DIR = PROJECT_DIR / "input"
OUT_DIR = PROJECT_DIR / "output"

CLINICAL_PATH = INPUT_DIR / "clinical_model3_paired72.tsv"
BASELINE_PATH = INPUT_DIR / "baseline_AAL3_12edges_KLDS_paired72.tsv"
EDGE_DEF_PATH = INPUT_DIR / "AAL3_Model3_12edge_definitions.tsv"

SEEDS = ["hippocampus", "amygdala", "dmpfc"]
TARGETS = ["ofc_medial_r", "ofc_anterior_r", "ofc_posterior_r", "insula_r"]
ALPHA = 0.05

FORMULA = (
    "MOCA_2 ~ MOCA_1 + UCLA_LS_1 + delta_UCLA_LS + KLDS_1 + "
    "delta_UCLA_LS*KLDS_1 + age_1 + sex + education_1 + FollowupInterval"
)

PREDICTORS = [
    "moca1_c",
    "ucla1_c",
    "delta_ucla_c",
    "baseline_klds_c",
    "interaction",
    "age1_c",
    "sex_1",
    "education1_c",
    "interval_c",
]


def centered(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    return values - values.mean()


def finite(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return np.nan
    return number if math.isfinite(number) else np.nan


def load_analysis_data() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    definitions = pd.read_csv(EDGE_DEF_PATH, sep="\t")
    definitions = definitions.loc[
        definitions["seed"].isin(SEEDS) & definitions["target"].isin(TARGETS)
    ].copy()
    definitions["seed_order"] = definitions["seed"].map(
        {value: index for index, value in enumerate(SEEDS)}
    )
    definitions["target_order"] = definitions["target"].map(
        {value: index for index, value in enumerate(TARGETS)}
    )
    definitions = definitions.sort_values(["seed_order", "target_order"]).drop(
        columns=["seed_order", "target_order"]
    )
    if len(definitions) != 12 or definitions["edge_key"].nunique() != 12:
        raise ValueError("Expected exactly 12 unique focused AAL3 edges")

    clinical_input = pd.read_csv(CLINICAL_PATH, sep="\t")
    clinical_columns = [
        "subject_id",
        "year_1",
        "year_2",
        "FollowupInterval",
        "age_1",
        "sex_1",
        "education_1",
        "UCLA_LS_1",
        "UCLA_LS_2",
        "delta_UCLA_LS",
        "MOCA_1",
        "MOCA_2",
    ]
    clinical = (
        clinical_input[clinical_columns]
        .drop_duplicates()
        .sort_values("subject_id")
        .reset_index(drop=True)
    )
    if len(clinical) != 72 or clinical["subject_id"].nunique() != 72:
        raise ValueError(f"Expected 72 unique clinical subjects, found {len(clinical)}")

    baseline = pd.read_csv(BASELINE_PATH, sep="\t").set_index("subject_id")

    clinical_ids = set(clinical["subject_id"])
    baseline_ids = set(baseline.index)
    if clinical_ids != baseline_ids:
        raise ValueError(
            "Clinical/baseline-AAL3 subject mismatch: "
            f"clinical-only={sorted(clinical_ids - baseline_ids)}, "
            f"AAL3-only={sorted(baseline_ids - clinical_ids)}"
        )

    frames: list[pd.DataFrame] = []
    for edge in definitions.to_dict(orient="records"):
        value_column = f"{edge['edge_key']}__KLDS"
        if value_column not in baseline:
            raise KeyError(f"Missing KLDS input column: {value_column}")
        frame = clinical.copy()
        frame["edge_key"] = edge["edge_key"]
        frame["edge"] = edge["edge"]
        frame["seed"] = edge["seed"]
        frame["target_region"] = edge["target"]
        frame["KLDS_1"] = frame["subject_id"].map(baseline[value_column])
        frames.append(frame)

    analysis = pd.concat(frames, ignore_index=True)
    required = [
        "MOCA_2",
        "MOCA_1",
        "UCLA_LS_1",
        "delta_UCLA_LS",
        "KLDS_1",
        "age_1",
        "sex_1",
        "education_1",
        "FollowupInterval",
    ]
    missing = analysis[required].isna().sum()
    if int(missing.sum()) != 0:
        raise ValueError(f"Missing Model 3 inputs:\n{missing}")

    pairing_qc = {
        "clinical_subjects": len(clinical_ids),
        "baseline_AAL3_subjects": int(len(baseline)),
        "analysis_subjects": len(baseline_ids),
        "focused_edges": int(definitions["edge_key"].nunique()),
        "analysis_rows": int(len(analysis)),
        "missing_model3_values": int(analysis[required].isna().sum().sum()),
        "followup_interval_counts": {
            str(int(key)): int(value)
            for key, value in clinical["FollowupInterval"]
            .value_counts()
            .sort_index()
            .items()
        },
    }
    return analysis, definitions, pairing_qc


def prepare_edge_frame(frame: pd.DataFrame) -> pd.DataFrame:
    work = frame.copy()
    work["moca1_c"] = centered(work["MOCA_1"])
    work["ucla1_c"] = centered(work["UCLA_LS_1"])
    work["delta_ucla_c"] = centered(work["delta_UCLA_LS"])
    work["baseline_klds_c"] = centered(work["KLDS_1"])
    work["interaction"] = work["delta_ucla_c"] * work["baseline_klds_c"]
    work["age1_c"] = centered(work["age_1"])
    work["education1_c"] = centered(work["education_1"])
    work["interval_c"] = centered(work["FollowupInterval"])
    return work


def fit_edge(
    frame: pd.DataFrame, edge: dict[str, object]
) -> tuple[
    dict[str, object],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    data = prepare_edge_frame(frame)
    needed = ["MOCA_2"] + PREDICTORS
    data = data[needed + ["delta_UCLA_LS", "KLDS_1"]].dropna().copy()
    y = data["MOCA_2"].astype(float)
    x = sm.add_constant(data[PREDICTORS].astype(float), has_constant="add")
    model = sm.OLS(y, x).fit()
    robust = model.get_robustcov_results(cov_type="HC3")

    target = "interaction"
    target_index = list(x.columns).index(target)
    beta = finite(model.params[target])
    se = finite(model.bse[target])
    t_value = finite(model.tvalues[target])
    p_value = finite(model.pvalues[target])
    ci = model.conf_int(alpha=0.05).loc[target]
    hc3_ci = robust.conf_int(alpha=0.05)[target_index]
    partial_r = np.sign(t_value) * math.sqrt(
        t_value**2 / (t_value**2 + model.df_resid)
    )
    standardized_beta = beta * data[target].std(ddof=1) / y.std(ddof=1)

    reduced_predictors = [name for name in PREDICTORS if name != target]
    reduced_x = sm.add_constant(
        data[reduced_predictors].astype(float), has_constant="add"
    )
    reduced_model = sm.OLS(y, reduced_x).fit()
    interaction_f, interaction_f_p, interaction_df = model.compare_f_test(reduced_model)

    vif_values = {
        name: finite(variance_inflation_factor(x.to_numpy(), index))
        for index, name in enumerate(x.columns)
        if name != "const"
    }
    residuals = np.asarray(model.resid)
    cooks = np.asarray(model.get_influence().cooks_distance[0])
    cook_threshold = 4.0 / len(data)

    target_row = {
        "edge_key": edge["edge_key"],
        "edge": edge["edge"],
        "seed": edge["seed"],
        "target_region": edge["target"],
        "n": int(model.nobs),
        "residual_df": int(model.df_resid),
        "interaction_beta": beta,
        "interaction_standardized_beta": finite(standardized_beta),
        "interaction_se": se,
        "interaction_ci95_low": finite(ci.iloc[0]),
        "interaction_ci95_high": finite(ci.iloc[1]),
        "interaction_t": t_value,
        "interaction_partial_r": finite(partial_r),
        "interaction_p_uncorrected": p_value,
        "interaction_hc3_se": finite(robust.bse[target_index]),
        "interaction_hc3_ci95_low": finite(hc3_ci[0]),
        "interaction_hc3_ci95_high": finite(hc3_ci[1]),
        "interaction_hc3_t": finite(robust.tvalues[target_index]),
        "interaction_hc3_p": finite(robust.pvalues[target_index]),
        "interaction_delta_r_squared": finite(model.rsquared - reduced_model.rsquared),
        "interaction_partial_f": finite(interaction_f),
        "interaction_partial_f_df": int(interaction_df),
        "interaction_partial_f_p": finite(interaction_f_p),
        "r_squared": finite(model.rsquared),
        "adjusted_r_squared": finite(model.rsquared_adj),
        "model_f": finite(model.fvalue),
        "model_f_p": finite(model.f_pvalue),
        "interaction_vif": vif_values[target],
        "max_predictor_vif": max(vif_values.values()),
        "shapiro_residual_p": finite(stats.shapiro(residuals).pvalue),
        "breusch_pagan_p": finite(het_breuschpagan(residuals, x)[1]),
        "max_cooks_d": finite(cooks.max()),
        "n_cooks_d_gt_4_over_n": int((cooks > cook_threshold).sum()),
        "formula": FORMULA,
    }

    coefficient_rows: list[dict[str, object]] = []
    conf_int = model.conf_int(alpha=0.05)
    robust_conf_int = robust.conf_int(alpha=0.05)
    y_sd = y.std(ddof=1)
    for index, term in enumerate(x.columns):
        x_sd = data[term].std(ddof=1) if term != "const" else np.nan
        standardized = (
            model.params[term] * x_sd / y_sd if term != "const" else np.nan
        )
        coefficient_rows.append(
            {
                "edge_key": edge["edge_key"],
                "seed": edge["seed"],
                "target_region": edge["target"],
                "term": term,
                "beta": finite(model.params[term]),
                "standardized_beta": finite(standardized),
                "se": finite(model.bse[term]),
                "ci95_low": finite(conf_int.loc[term, 0]),
                "ci95_high": finite(conf_int.loc[term, 1]),
                "t": finite(model.tvalues[term]),
                "p_uncorrected": finite(model.pvalues[term]),
                "hc3_se": finite(robust.bse[index]),
                "hc3_ci95_low": finite(robust_conf_int[index, 0]),
                "hc3_ci95_high": finite(robust_conf_int[index, 1]),
                "hc3_t": finite(robust.tvalues[index]),
                "hc3_p": finite(robust.pvalues[index]),
            }
        )

    slopes: list[dict[str, object]] = []
    klds_sd = frame["KLDS_1"].std(ddof=1)
    covariance = model.cov_params()
    critical = stats.t.ppf(0.975, model.df_resid)
    for label, klds_centered in [
        ("Low_-1SD", -klds_sd),
        ("Mean", 0.0),
        ("High_+1SD", klds_sd),
    ]:
        slope = (
            model.params["delta_ucla_c"]
            + klds_centered * model.params["interaction"]
        )
        variance = (
            covariance.loc["delta_ucla_c", "delta_ucla_c"]
            + klds_centered**2 * covariance.loc["interaction", "interaction"]
            + 2
            * klds_centered
            * covariance.loc["delta_ucla_c", "interaction"]
        )
        slope_se = math.sqrt(variance)
        slope_t = slope / slope_se
        slope_p = 2 * stats.t.sf(abs(slope_t), model.df_resid)
        slopes.append(
            {
                "edge_key": edge["edge_key"],
                "seed": edge["seed"],
                "target_region": edge["target"],
                "baseline_KLDS_level": label,
                "centered_KLDS_value": klds_centered,
                "delta_UCLA_slope_on_MOCA2": finite(slope),
                "se": finite(slope_se),
                "ci95_low": finite(slope - critical * slope_se),
                "ci95_high": finite(slope + critical * slope_se),
                "t": finite(slope_t),
                "p_uncorrected": finite(slope_p),
            }
        )

    return target_row, coefficient_rows, slopes


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    analysis, definitions, pairing_qc = load_analysis_data()

    target_rows: list[dict[str, object]] = []
    coefficient_rows: list[dict[str, object]] = []
    slope_rows: list[dict[str, object]] = []
    for edge in definitions.to_dict(orient="records"):
        frame = analysis.loc[analysis["edge_key"].eq(edge["edge_key"])].copy()
        if len(frame) != 72 or frame["subject_id"].nunique() != 72:
            raise ValueError(f"{edge['edge_key']}: expected 72 unique participants")
        target_row, coefficients, slopes = fit_edge(frame, edge)
        target_rows.append(target_row)
        coefficient_rows.extend(coefficients)
        slope_rows.extend(slopes)

    results = pd.DataFrame(target_rows)
    reject, q_values, _, _ = multipletests(
        results["interaction_p_uncorrected"].to_numpy(),
        alpha=ALPHA,
        method="fdr_bh",
    )
    hc3_reject, hc3_q_values, _, _ = multipletests(
        results["interaction_hc3_p"].to_numpy(),
        alpha=ALPHA,
        method="fdr_bh",
    )
    results["interaction_q_BH_FDR_12"] = q_values
    results["interaction_hc3_q_BH_FDR_12"] = hc3_q_values
    results["interaction_significant_uncorrected"] = (
        results["interaction_p_uncorrected"] < ALPHA
    )
    results["interaction_significant_BH_FDR_12"] = reject
    results["interaction_hc3_significant_BH_FDR_12"] = hc3_reject

    edge_order = {
        key: index for index, key in enumerate(definitions["edge_key"].tolist())
    }
    results["_edge_order"] = results["edge_key"].map(edge_order)
    results = results.sort_values("_edge_order").drop(columns="_edge_order")

    coefficients = pd.DataFrame(coefficient_rows)
    coefficients["_edge_order"] = coefficients["edge_key"].map(edge_order)
    coefficients = coefficients.sort_values(["_edge_order", "term"]).drop(
        columns="_edge_order"
    )
    slopes = pd.DataFrame(slope_rows)
    slope_order = {"Low_-1SD": 0, "Mean": 1, "High_+1SD": 2}
    slopes["_edge_order"] = slopes["edge_key"].map(edge_order)
    slopes["_slope_order"] = slopes["baseline_KLDS_level"].map(slope_order)
    slopes = slopes.sort_values(["_edge_order", "_slope_order"]).drop(
        columns=["_edge_order", "_slope_order"]
    )

    output_columns = [
        "subject_id",
        "edge_key",
        "edge",
        "seed",
        "target_region",
        "KLDS_1",
        "UCLA_LS_1",
        "UCLA_LS_2",
        "delta_UCLA_LS",
        "MOCA_1",
        "MOCA_2",
        "age_1",
        "sex_1",
        "education_1",
        "FollowupInterval",
    ]
    analysis_out = analysis[output_columns].copy()
    analysis_out["_edge_order"] = analysis_out["edge_key"].map(edge_order)
    analysis_out = analysis_out.sort_values(["_edge_order", "subject_id"]).drop(
        columns="_edge_order"
    )

    qc = {
        "analysis": "AAL3 KLDS Model 3-only focused analysis",
        "formula": FORMULA,
        "target_term": "centered delta_UCLA_LS x centered baseline KLDS",
        "paired_subjects": 72,
        "seeds": SEEDS,
        "targets": TARGETS,
        "edges": 12,
        "primary_tests": 12,
        "primary_correction": "BH-FDR across the 12 interaction terms",
        "common_covariates": [
            "baseline MOCA",
            "baseline UCLA-LS",
            "baseline age",
            "sex",
            "baseline education",
            "FollowupInterval",
        ],
        "lower_order_terms": ["delta_UCLA_LS", "baseline KLDS"],
        "primary_standard_errors": "conventional OLS",
        "sensitivity_standard_errors": "HC3",
        "pairing": pairing_qc,
        "fdr_significant_interactions": int(
            results["interaction_significant_BH_FDR_12"].sum()
        ),
        "raw_p_significant_interactions": int(
            results["interaction_significant_uncorrected"].sum()
        ),
        "hc3_fdr_significant_interactions": int(
            results["interaction_hc3_significant_BH_FDR_12"].sum()
        ),
        "missing_output_values": int(
            results[
                [
                    "interaction_beta",
                    "interaction_se",
                    "interaction_t",
                    "interaction_p_uncorrected",
                    "interaction_q_BH_FDR_12",
                ]
            ]
            .isna()
            .sum()
            .sum()
        ),
    }

    results.to_csv(
        OUT_DIR / "AAL3_KLDS_Model3_12edges_interaction_results.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    coefficients.to_csv(
        OUT_DIR / "AAL3_KLDS_Model3_12edges_all_coefficients.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    slopes.to_csv(
        OUT_DIR / "AAL3_KLDS_Model3_12edges_simple_slopes.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    analysis_out.to_csv(
        OUT_DIR / "AAL3_KLDS_Model3_12edges_paired72_analysis_data.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    definitions.to_csv(
        OUT_DIR / "AAL3_KLDS_Model3_12edge_definitions.tsv",
        sep="\t",
        index=False,
    )
    with (OUT_DIR / "AAL3_KLDS_Model3_12edges_QC.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)

    method = f"""AAL3 KLDS MODEL 3-ONLY ANALYSIS

Sample
- 72 participants with complete baseline KLDS and baseline/follow-up clinical data.
- FollowupInterval is included because observed intervals are 2, 3, or 4 years.

Edges
- Three bilateral-mean seeds: hippocampus, amygdala, dmPFC.
- Four right targets: medial OFC, anterior OFC, posterior OFC, insula.
- Twelve prespecified KLDS edges.

Model
{FORMULA}

- Continuous predictors are mean-centered before forming the interaction.
- Sex retains the source 0/1 coding.
- Target: centered delta_UCLA_LS x centered baseline KLDS interaction.

Inference
- Ordinary least squares, conventional two-sided p-values as primary inference.
- BH-FDR across the 12 interaction tests.
- HC3 standard errors, p-values, and a separate 12-test BH-FDR are reported as
  sensitivity statistics.
- Simple slopes describe delta_UCLA_LS at baseline KLDS mean -1 SD, mean, and mean +1 SD.
"""
    (OUT_DIR / "AAL3_KLDS_Model3_12edges_METHOD.txt").write_text(
        method, encoding="utf-8"
    )

    print(json.dumps(qc, ensure_ascii=False, indent=2))
    print("\nMODEL 3 INTERACTION RESULTS")
    print(
        results[
            [
                "edge_key",
                "interaction_partial_r",
                "interaction_beta",
                "interaction_se",
                "interaction_ci95_low",
                "interaction_ci95_high",
                "interaction_t",
                "interaction_p_uncorrected",
                "interaction_q_BH_FDR_12",
                "interaction_hc3_p",
                "interaction_hc3_q_BH_FDR_12",
                "interaction_delta_r_squared",
                "adjusted_r_squared",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
