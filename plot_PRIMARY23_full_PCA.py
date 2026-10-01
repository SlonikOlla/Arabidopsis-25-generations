#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
PCADIR = ROOT / "results/wgs_variants/PRIMARY23_PCA"
FIGDIR = PCADIR / "figures"

META = ROOT / "metadata/wgs_primary23_groups.tsv"

ALL_PREFIX = PCADIR / "PRIMARY23_all23_5kb_full"
NO65_PREFIX = PCADIR / "PRIMARY23_without_SRR17867665_5kb_full"

ALL_ANNOTATED = PCADIR / "PRIMARY23_all23_5kb_full_annotated.tsv"
NO65_ANNOTATED = (
    PCADIR /
    "PRIMARY23_without_SRR17867665_5kb_full_annotated.tsv"
)

ALL_SUMMARY = PCADIR / "PRIMARY23_all23_5kb_full_group_summary.tsv"
NO65_SUMMARY = (
    PCADIR /
    "PRIMARY23_without_SRR17867665_5kb_full_group_summary.tsv"
)

FIG_PNG = FIGDIR / "PRIMARY23_PCA_primary_and_sensitivity.png"
FIG_PDF = FIGDIR / "PRIMARY23_PCA_primary_and_sensitivity.pdf"

GROUP_ORDER = ["F2C", "F25C", "F25H", "F25Cd", "F25UV"]

COLORS = {
    "F2C":   "#7F7F7F",
    "F25C":  "#009E73",
    "F25H":  "#D55E00",
    "F25Cd": "#0072B2",
    "F25UV": "#CC79A7",
}

MARKERS = {
    "F2C":   "o",
    "F25C":  "s",
    "F25H":  "^",
    "F25Cd": "D",
    "F25UV": "P",
}


def read_variance(prefix):
    table = pd.read_csv(
        str(prefix) + ".variance.tsv",
        sep="\t"
    )

    required = {"component", "eigenvalue", "variance_percent"}

    if not required.issubset(table.columns):
        raise RuntimeError(
            f"Malformed variance table for {prefix}"
        )

    return dict(
        zip(table["component"], table["variance_percent"])
    )


def read_and_annotate(prefix, expected_samples, excluded=None):
    eigenvec = pd.read_csv(
        str(prefix) + ".eigenvec",
        sep=r"\s+"
    )

    eigenvec.columns = [
        column.lstrip("#") for column in eigenvec.columns
    ]

    required = {"FID", "IID", "PC1", "PC2", "PC3"}

    if not required.issubset(eigenvec.columns):
        raise RuntimeError(
            f"Required PCA columns absent from {prefix}.eigenvec"
        )

    metadata = pd.read_csv(META, sep=r"\s+")

    if set(metadata.columns) != {"sample", "group"}:
        raise RuntimeError("Malformed group metadata table")

    if metadata["sample"].duplicated().any():
        raise RuntimeError("Duplicate samples in group metadata")

    annotated = eigenvec.merge(
        metadata,
        left_on="IID",
        right_on="sample",
        how="left",
        validate="one_to_one"
    )

    if annotated["group"].isna().any():
        missing = annotated.loc[
            annotated["group"].isna(), "IID"
        ].tolist()
        raise RuntimeError(
            f"Samples without group annotation: {missing}"
        )

    if len(annotated) != expected_samples:
        raise RuntimeError(
            f"{prefix}: expected {expected_samples} samples; "
            f"found {len(annotated)}"
        )

    if annotated["IID"].nunique() != expected_samples:
        raise RuntimeError(
            f"{prefix}: duplicated PCA sample identifiers"
        )

    observed_groups = set(annotated["group"])
    expected_groups = set(GROUP_ORDER)

    if observed_groups != expected_groups:
        raise RuntimeError(
            f"{prefix}: group mismatch: {observed_groups}"
        )

    if excluded is not None and excluded in set(annotated["IID"]):
        raise RuntimeError(
            f"{excluded} remains in sensitivity PCA"
        )

    annotated["group"] = pd.Categorical(
        annotated["group"],
        categories=GROUP_ORDER,
        ordered=True
    )

    annotated = annotated.sort_values(
        ["group", "IID"]
    ).reset_index(drop=True)

    return annotated


def write_group_summary(data, output):
    rows = []

    for group in GROUP_ORDER:
        subset = data.loc[data["group"] == group]

        row = {
            "group": group,
            "n": len(subset),
        }

        for pc in ["PC1", "PC2", "PC3"]:
            row[f"{pc}_mean"] = subset[pc].mean()
            row[f"{pc}_sd"] = subset[pc].std(ddof=1)

        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(output, sep="\t", index=False)
    return summary


def padded_limits(series):
    minimum = float(series.min())
    maximum = float(series.max())
    span = maximum - minimum

    if span == 0:
        span = 1.0

    padding = span * 0.12
    return minimum - padding, maximum + padding


def plot_panel(
    axis,
    data,
    x_pc,
    y_pc,
    variance,
    title,
    panel_letter,
    label_outlier=False
):
    for group in GROUP_ORDER:
        subset = data.loc[data["group"] == group]

        axis.scatter(
            subset[x_pc],
            subset[y_pc],
            s=85,
            marker=MARKERS[group],
            color=COLORS[group],
            edgecolor="black",
            linewidth=0.7,
            alpha=0.92,
            zorder=3
        )

    axis.axhline(
        0,
        color="#D0D0D0",
        linewidth=0.8,
        zorder=1
    )
    axis.axvline(
        0,
        color="#D0D0D0",
        linewidth=0.8,
        zorder=1
    )

    axis.set_xlabel(
        f"{x_pc} ({variance[x_pc]:.2f}%)",
        fontsize=11
    )
    axis.set_ylabel(
        f"{y_pc} ({variance[y_pc]:.2f}%)",
        fontsize=11
    )

    axis.set_title(title, fontsize=12, pad=10)
    axis.text(
        -0.12,
        1.06,
        panel_letter,
        transform=axis.transAxes,
        fontsize=15,
        fontweight="bold",
        va="top"
    )

    axis.grid(
        True,
        color="#ECECEC",
        linewidth=0.7,
        zorder=0
    )

    axis.set_xlim(padded_limits(data[x_pc]))
    axis.set_ylim(padded_limits(data[y_pc]))

    axis.tick_params(labelsize=9)

    for spine in axis.spines.values():
        spine.set_linewidth(0.8)

    if label_outlier:
        outlier = data.loc[
            data["IID"] == "SRR17867665"
        ]

        if len(outlier) == 1:
            x_value = float(outlier.iloc[0][x_pc])
            y_value = float(outlier.iloc[0][y_pc])

            axis.annotate(
                "SRR17867665",
                xy=(x_value, y_value),
                xytext=(8, 8),
                textcoords="offset points",
                fontsize=8,
                color=COLORS["F25H"],
                arrowprops={
                    "arrowstyle": "-",
                    "color": COLORS["F25H"],
                    "linewidth": 0.8,
                },
                zorder=4
            )


FIGDIR.mkdir(parents=True, exist_ok=True)

all_data = read_and_annotate(
    ALL_PREFIX,
    expected_samples=23
)

no65_data = read_and_annotate(
    NO65_PREFIX,
    expected_samples=22,
    excluded="SRR17867665"
)

all_variance = read_variance(ALL_PREFIX)
no65_variance = read_variance(NO65_PREFIX)

all_data.to_csv(
    ALL_ANNOTATED,
    sep="\t",
    index=False
)

no65_data.to_csv(
    NO65_ANNOTATED,
    sep="\t",
    index=False
)

all_summary = write_group_summary(
    all_data,
    ALL_SUMMARY
)

no65_summary = write_group_summary(
    no65_data,
    NO65_SUMMARY
)

print("===== NUMERICAL AUDIT =====")
print(f"Primary PCA samples:     {len(all_data)}")
print(f"Sensitivity PCA samples: {len(no65_data)}")
print(f"Primary unique samples:  {all_data['IID'].nunique()}")
print(f"Sensitivity unique:      {no65_data['IID'].nunique()}")
print(
    "SRR17867665 in primary:",
    "SRR17867665" in set(all_data["IID"])
)
print(
    "SRR17867665 in sensitivity:",
    "SRR17867665" in set(no65_data["IID"])
)

print()
print("===== GROUP COUNTS: PRIMARY =====")
print(
    all_data["group"]
    .value_counts(sort=False)
    .to_string()
)

print()
print("===== GROUP COUNTS: SENSITIVITY =====")
print(
    no65_data["group"]
    .value_counts(sort=False)
    .to_string()
)

for label, data in [
    ("PRIMARY", all_data),
    ("SENSITIVITY", no65_data),
]:
    print()
    print(f"===== COORDINATE RANGES: {label} =====")

    for pc in ["PC1", "PC2", "PC3"]:
        minimum_row = data.loc[data[pc].idxmin()]
        maximum_row = data.loc[data[pc].idxmax()]

        print(
            f"{pc}: "
            f"min={minimum_row[pc]:.8f} "
            f"({minimum_row['IID']}); "
            f"max={maximum_row[pc]:.8f} "
            f"({maximum_row['IID']})"
        )

fig, axes = plt.subplots(
    2,
    2,
    figsize=(12.0, 10.0)
)

plot_panel(
    axes[0, 0],
    all_data,
    "PC1",
    "PC2",
    all_variance,
    "All 23 samples",
    "A",
    label_outlier=True
)

plot_panel(
    axes[0, 1],
    no65_data,
    "PC1",
    "PC2",
    no65_variance,
    "Sensitivity analysis without SRR17867665",
    "B"
)

plot_panel(
    axes[1, 0],
    all_data,
    "PC1",
    "PC3",
    all_variance,
    "All 23 samples",
    "C",
    label_outlier=True
)

plot_panel(
    axes[1, 1],
    no65_data,
    "PC1",
    "PC3",
    no65_variance,
    "Sensitivity analysis without SRR17867665",
    "D"
)

legend_handles = [
    Line2D(
        [0],
        [0],
        marker=MARKERS[group],
        linestyle="none",
        markerfacecolor=COLORS[group],
        markeredgecolor="black",
        markeredgewidth=0.7,
        markersize=9,
        label=group
    )
    for group in GROUP_ORDER
]

fig.legend(
    handles=legend_handles,
    labels=GROUP_ORDER,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.985),
    ncol=5,
    frameon=False,
    fontsize=10
)

fig.suptitle(
    "Genetic structure of 25-generation Arabidopsis stress lineages",
    fontsize=15,
    y=0.995
)

fig.subplots_adjust(
    top=0.91,
    bottom=0.08,
    left=0.09,
    right=0.98,
    hspace=0.28,
    wspace=0.24
)

fig.savefig(
    FIG_PNG,
    dpi=400,
    bbox_inches="tight"
)

fig.savefig(
    FIG_PDF,
    bbox_inches="tight"
)

plt.close(fig)

for output in [
    ALL_ANNOTATED,
    NO65_ANNOTATED,
    ALL_SUMMARY,
    NO65_SUMMARY,
    FIG_PNG,
    FIG_PDF,
]:
    if not output.exists() or output.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty output: {output}")

print()
print("===== GROUP SUMMARY: PRIMARY =====")
print(all_summary.to_string(index=False))

print()
print("===== GROUP SUMMARY: SENSITIVITY =====")
print(no65_summary.to_string(index=False))

print()
print("===== OUTPUTS =====")
print(ALL_ANNOTATED)
print(NO65_ANNOTATED)
print(ALL_SUMMARY)
print(NO65_SUMMARY)
print(FIG_PNG)
print(FIG_PDF)
print("PRIMARY23 PCA TABLES AND FIGURE: COMPLETE")
