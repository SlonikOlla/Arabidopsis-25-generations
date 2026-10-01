#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as t_distribution

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
BASEDIR = ROOT / "results/wgbs_PRIMARY22_500bp"
OUTDIR = BASEDIR / "UVC_background_adjusted"
DESIGN_FILE = (
    ROOT /
    "metadata/wgbs_F25C_F25UV_background_adjusted_design.tsv"
)

CONTEXTS = ["CG", "CHG", "CHH"]

EFFECT_THRESHOLDS = {
    "CG": 0.075,
    "CHG": 0.075,
    "CHH": 0.050,
}

MATRIX_FILES = {
    context:
    BASEDIR /
    f"PRIMARY22_{context}_500bp_weighted_methylation.tsv.gz"
    for context in CONTEXTS
}

OUTDIR.mkdir(parents=True, exist_ok=True)

design = pd.read_csv(
    DESIGN_FILE,
    sep="\t",
    dtype=str,
)

expected_samples = 10

if len(design) != expected_samples:
    raise RuntimeError(
        f"Expected 10 design samples, found {len(design)}"
    )

if design["sample"].duplicated().any():
    raise RuntimeError("Duplicate samples in design")

design["treatment_code"] = (
    design["treatment"] == "UVC"
).astype(float)

design["background_code"] = (
    design["background"] == "B"
).astype(float)

sample_order = design["sample"].tolist()

X = np.column_stack([
    np.ones(expected_samples),
    design["treatment_code"].to_numpy(),
    design["background_code"].to_numpy(),
])

rank = np.linalg.matrix_rank(X)

if rank != 3:
    raise RuntimeError(f"Design matrix rank is {rank}, expected 3")

XtX_inverse = np.linalg.inv(X.T @ X)
coefficient_variance_factor = XtX_inverse[1, 1]

summary_rows = []

for context in CONTEXTS:

    matrix_file = MATRIX_FILES[context]

    print(f"===== {context} =====", flush=True)
    print(f"READING: {matrix_file}", flush=True)

    table = pd.read_csv(
        matrix_file,
        sep="\t",
        dtype={"chrom": str},
    )

    required = {
        "locus",
        *sample_order,
    }

    missing_columns = required.difference(table.columns)

    if missing_columns:
        raise RuntimeError(
            f"{context}: missing columns: "
            + ",".join(sorted(missing_columns))
        )

    complete = table[sample_order].notna().all(axis=1)
    tested = table.loc[complete].copy()

    Y = tested[sample_order].to_numpy(
        dtype=float,
        copy=True,
    ).T

    beta = XtX_inverse @ X.T @ Y

    fitted = X @ beta
    residuals = Y - fitted

    degrees_freedom = expected_samples - rank

    residual_sum_squares = np.sum(
        residuals ** 2,
        axis=0,
    )

    residual_variance = (
        residual_sum_squares /
        degrees_freedom
    )

    treatment_standard_error = np.sqrt(
        residual_variance *
        coefficient_variance_factor
    )

    treatment_beta = beta[1, :]

    with np.errstate(divide="ignore", invalid="ignore"):
        t_statistic = (
            treatment_beta /
            treatment_standard_error
        )

    p_value = 2.0 * t_distribution.sf(
        np.abs(t_statistic),
        df=degrees_freedom,
    )

    zero_variance = treatment_standard_error == 0

    p_value[
        zero_variance &
        (treatment_beta == 0)
    ] = 1.0

    p_value[
        zero_variance &
        (treatment_beta != 0)
    ] = 0.0

    p_value = np.where(
        np.isfinite(p_value),
        p_value,
        1.0,
    )

    order = np.argsort(p_value)
    ranked = p_value[order]
    number_tests = len(ranked)

    adjusted_ranked = (
        ranked *
        number_tests /
        np.arange(1, number_tests + 1)
    )

    adjusted_ranked = np.minimum.accumulate(
        adjusted_ranked[::-1]
    )[::-1]

    adjusted_ranked = np.minimum(
        adjusted_ranked,
        1.0,
    )

    fdr = np.empty(number_tests, dtype=float)
    fdr[order] = adjusted_ranked

    effect_threshold = EFFECT_THRESHOLDS[context]

    direction = np.full(
        number_tests,
        "NONE",
        dtype=object,
    )

    significant = (
        (fdr < 0.05) &
        (np.abs(treatment_beta) >= effect_threshold)
    )

    direction[
        significant &
        (treatment_beta > 0)
    ] = "UVC_HYPER"

    direction[
        significant &
        (treatment_beta < 0)
    ] = "UVC_HYPO"

    control_a = design.loc[
        (design["treatment"] == "CONTROL") &
        (design["background"] == "A"),
        "sample"
    ].tolist()

    uvc_a = design.loc[
        (design["treatment"] == "UVC") &
        (design["background"] == "A"),
        "sample"
    ].tolist()

    control_b = design.loc[
        (design["treatment"] == "CONTROL") &
        (design["background"] == "B"),
        "sample"
    ].tolist()

    uvc_b = design.loc[
        (design["treatment"] == "UVC") &
        (design["background"] == "B"),
        "sample"
    ].tolist()

    coordinates = tested["locus"].str.extract(
        r"^(?P<chrom>[^:]+):"
        r"(?P<start>\d+)-(?P<end>\d+)$"
    )

    if coordinates.isna().any().any():
        bad_loci = tested.loc[
            coordinates.isna().any(axis=1),
            "locus"
        ].head(10).tolist()

        raise RuntimeError(
            f"{context}: malformed loci: {bad_loci}"
        )

    results = pd.DataFrame({
        "chrom": coordinates["chrom"].to_numpy(),
        "start": coordinates["start"].astype(int).to_numpy(),
        "end": coordinates["end"].astype(int).to_numpy(),
    })

    results["context"] = context
    results["n_samples"] = expected_samples

    results["CONTROL_A_mean"] = (
        tested[control_a].mean(axis=1).to_numpy()
    )

    results["UVC_A_mean"] = (
        tested[uvc_a].mean(axis=1).to_numpy()
    )

    results["delta_A"] = (
        results["UVC_A_mean"] -
        results["CONTROL_A_mean"]
    )

    results["CONTROL_B_mean"] = (
        tested[control_b].mean(axis=1).to_numpy()
    )

    results["UVC_B_mean"] = (
        tested[uvc_b].mean(axis=1).to_numpy()
    )

    results["delta_B"] = (
        results["UVC_B_mean"] -
        results["CONTROL_B_mean"]
    )

    results["UVC_treatment_beta"] = treatment_beta
    results["background_B_beta"] = beta[2, :]
    results["treatment_SE"] = treatment_standard_error
    results["t_statistic"] = t_statistic
    results["P_value"] = p_value
    results["FDR"] = fdr
    results["effect_threshold"] = effect_threshold
    results["direction"] = direction

    results = results.sort_values(
        ["FDR", "P_value", "chrom", "start"]
    )

    full_file = (
        OUTDIR /
        f"PRIMARY22_{context}_UVC_background_adjusted_all.tsv.gz"
    )

    significant_file = (
        OUTDIR /
        f"PRIMARY22_{context}_UVC_background_adjusted_significant.tsv"
    )

    results.to_csv(
        full_file,
        sep="\t",
        index=False,
        compression="gzip",
    )

    results.loc[
        results["direction"] != "NONE"
    ].to_csv(
        significant_file,
        sep="\t",
        index=False,
    )

    hyper = int(
        np.sum(direction == "UVC_HYPER")
    )

    hypo = int(
        np.sum(direction == "UVC_HYPO")
    )

    fdr_only = int(
        np.sum(fdr < 0.05)
    )

    concordant_direction = int(
        np.sum(
            ((results["delta_A"] > 0) &
             (results["delta_B"] > 0)) |
            ((results["delta_A"] < 0) &
             (results["delta_B"] < 0))
        )
    )

    summary_rows.append({
        "context": context,
        "input_windows": len(table),
        "complete_10_sample_windows": len(results),
        "degrees_freedom": degrees_freedom,
        "effect_threshold": effect_threshold,
        "FDR_below_0.05": fdr_only,
        "significant_hypermethylated": hyper,
        "significant_hypomethylated": hypo,
        "significant_total": hyper + hypo,
        "windows_with_same_delta_direction":
            concordant_direction,
    })

    print(f"INPUT WINDOWS: {len(table)}")
    print(f"COMPLETE WINDOWS: {len(results)}")
    print(f"FDR < 0.05: {fdr_only}")
    print(f"UVC HYPER: {hyper}")
    print(f"UVC HYPO: {hypo}")

summary = pd.DataFrame(summary_rows)

summary_file = (
    OUTDIR /
    "PRIMARY22_UVC_background_adjusted_summary.tsv"
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False,
)

print("===== FINAL SUMMARY =====")
print(summary.to_string(index=False))
print("PRIMARY22 UVC BACKGROUND-ADJUSTED TEST: COMPLETE")
