#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as t_distribution

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")

BASEDIR = (
    ROOT /
    "results/wgbs_PRIMARY22_500bp"
)

OUTDIR = (
    BASEDIR /
    "UVC_background_adjusted"
)

DESIGN_FILE = (
    ROOT /
    "metadata/wgbs_F25C_F25UV_background_adjusted_design.tsv"
)

CONTEXTS = ["CG", "CHH"]

design = pd.read_csv(
    DESIGN_FILE,
    sep="\t",
    dtype=str,
)

design["T"] = (
    design["treatment"] == "UVC"
).astype(float)

design["B"] = (
    design["background"] == "B"
).astype(float)

design["T_by_B"] = (
    design["T"] *
    design["B"]
)

samples = design["sample"].tolist()

X = np.column_stack([
    np.ones(len(design)),
    design["T"].to_numpy(),
    design["B"].to_numpy(),
    design["T_by_B"].to_numpy(),
])

rank = np.linalg.matrix_rank(X)

if rank != 4:
    raise RuntimeError(
        f"Interaction design rank={rank}, expected 4"
    )

inverse = np.linalg.inv(X.T @ X)
interaction_variance_factor = inverse[3, 3]
degrees_freedom = len(design) - rank

all_results = []

for context in CONTEXTS:

    print(f"===== {context} =====", flush=True)

    candidate_file = (
        OUTDIR /
        f"PRIMARY22_{context}_"
        "UVC_background_adjusted_significant.tsv"
    )

    matrix_file = (
        BASEDIR /
        f"PRIMARY22_{context}_"
        "500bp_weighted_methylation.tsv.gz"
    )

    candidates = pd.read_csv(
        candidate_file,
        sep="\t",
        dtype={"chrom": str},
    )

    matrix = pd.read_csv(
        matrix_file,
        sep="\t",
        dtype={"locus": str},
    ).set_index("locus")

    candidate_loci = (
        candidates["chrom"].astype(str) + ":" +
        candidates["start"].astype(str) + "-" +
        candidates["end"].astype(str)
    )

    missing = [
        locus for locus in candidate_loci
        if locus not in matrix.index
    ]

    if missing:
        raise RuntimeError(
            f"{context}: missing candidate loci: "
            + ",".join(missing[:10])
        )

    values = matrix.loc[
        candidate_loci,
        samples,
    ]

    if values.isna().any().any():
        raise RuntimeError(
            f"{context}: candidate matrix contains missing values"
        )

    Y = values.to_numpy(
        dtype=float,
        copy=True,
    ).T

    beta = inverse @ X.T @ Y
    fitted = X @ beta
    residuals = Y - fitted

    residual_sum_squares = np.sum(
        residuals ** 2,
        axis=0,
    )

    residual_variance = (
        residual_sum_squares /
        degrees_freedom
    )

    interaction_se = np.sqrt(
        residual_variance *
        interaction_variance_factor
    )

    interaction_beta = beta[3, :]

    with np.errstate(
        divide="ignore",
        invalid="ignore",
    ):
        interaction_t = (
            interaction_beta /
            interaction_se
        )

    interaction_p = (
        2.0 *
        t_distribution.sf(
            np.abs(interaction_t),
            df=degrees_freedom,
        )
    )

    zero_se = interaction_se == 0

    interaction_p[
        zero_se &
        (interaction_beta == 0)
    ] = 1.0

    interaction_p[
        zero_se &
        (interaction_beta != 0)
    ] = 0.0

    interaction_p = np.where(
        np.isfinite(interaction_p),
        interaction_p,
        1.0,
    )

    result = candidates.copy()
    result["locus"] = candidate_loci
    result["interaction_beta_deltaB_minus_deltaA"] = (
        interaction_beta
    )
    result["interaction_SE"] = interaction_se
    result["interaction_t"] = interaction_t
    result["interaction_P"] = interaction_p
    result["full_model_delta_A"] = beta[1, :]
    result["full_model_delta_B"] = (
        beta[1, :] +
        beta[3, :]
    )
    result["interaction_df"] = degrees_freedom

    result["same_direction"] = np.where(
        (
            (result["full_model_delta_A"] > 0) &
            (result["full_model_delta_B"] > 0)
        ) |
        (
            (result["full_model_delta_A"] < 0) &
            (result["full_model_delta_B"] < 0)
        ),
        "YES",
        "NO",
    )

    all_results.append(result)

combined = pd.concat(
    all_results,
    ignore_index=True,
)

def bh_adjust(values):

    values = np.asarray(
        values,
        dtype=float,
    )

    order = np.argsort(values)
    ranked = values[order]
    n = len(values)

    adjusted = (
        ranked *
        n /
        np.arange(1, n + 1)
    )

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.minimum(
        adjusted,
        1.0,
    )

    output = np.empty(n)
    output[order] = adjusted

    return output

combined["interaction_BH_FDR_all91"] = bh_adjust(
    combined["interaction_P"]
)

combined["interaction_BH_FDR_within_context"] = np.nan

for context in CONTEXTS:

    index = combined.index[
        combined["context"] == context
    ]

    combined.loc[
        index,
        "interaction_BH_FDR_within_context"
    ] = bh_adjust(
        combined.loc[index, "interaction_P"]
    )

combined["interaction_interpretation"] = "NO_EVIDENCE"

combined.loc[
    combined["interaction_P"] < 0.05,
    "interaction_interpretation"
] = "NOMINAL_BACKGROUND_HETEROGENEITY"

combined.loc[
    combined["interaction_BH_FDR_all91"] < 0.05,
    "interaction_interpretation"
] = "FDR_SUPPORTED_BACKGROUND_HETEROGENEITY"

outfile = (
    OUTDIR /
    "PRIMARY22_UVC_DMR_background_interaction_sensitivity.tsv"
)

combined.to_csv(
    outfile,
    sep="\t",
    index=False,
)

summary_rows = []

for context in CONTEXTS:

    subset = combined[
        combined["context"] == context
    ]

    summary_rows.append({
        "context": context,
        "DMRs": len(subset),
        "same_direction_A_and_B":
            int((subset["same_direction"] == "YES").sum()),
        "nominal_interaction_P_below_0.05":
            int((subset["interaction_P"] < 0.05).sum()),
        "interaction_FDR_below_0.05_all91":
            int((
                subset["interaction_BH_FDR_all91"] < 0.05
            ).sum()),
        "maximum_absolute_interaction":
            subset[
                "interaction_beta_deltaB_minus_deltaA"
            ].abs().max(),
        "median_absolute_interaction":
            subset[
                "interaction_beta_deltaB_minus_deltaA"
            ].abs().median(),
    })

summary = pd.DataFrame(summary_rows)

summary_file = (
    OUTDIR /
    "PRIMARY22_UVC_DMR_background_interaction_summary.tsv"
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False,
)

print("===== INTERACTION SUMMARY =====")
print(summary.to_string(index=False))
print("TOTAL DMRs:", len(combined))
print(
    "DIRECTION-DISCORDANT:",
    int((combined["same_direction"] == "NO").sum()),
)
print("WROTE:", outfile)
print("WROTE:", summary_file)
print("UVC DMR INTERACTION SENSITIVITY: COMPLETE")
