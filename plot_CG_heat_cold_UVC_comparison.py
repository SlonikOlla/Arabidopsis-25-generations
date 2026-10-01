#!/usr/bin/env python3

from pathlib import Path
import csv
import math

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
BASE = ROOT / "results/cross_stress_CG_comparison"
FIGDIR = BASE / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

MASTER = BASE / "CG_heat_cold_UVC_switch_master.tsv"
ENRICH = BASE / "CG_switch_overlap_enrichment.tsv"
RECURRENT = (
    BASE /
    "CG_recurrent_UVC_heat_cold_gene_descriptions.tsv"
)

PNG = FIGDIR / "CG_heat_cold_UVC_comparison.png"
PDF = FIGDIR / "CG_heat_cold_UVC_comparison.pdf"
AUDIT = BASE / "CG_heat_cold_UVC_figure_source_audit.tsv"

COLORS = {
    "GAIN": "#D55E00",
    "LOSS": "#0072B2",
    "ABSENT": "#E6E6E6",
}

STRESSES = ["HEAT", "COLD", "UVC"]


master = pd.read_csv(
    MASTER,
    sep="\t",
    dtype={"chrom": str}
)

enrich = pd.read_csv(
    ENRICH,
    sep="\t"
)

recurrent = pd.read_csv(
    RECURRENT,
    sep="\t",
    dtype={"chrom": str}
)


# ------------------------------------------------------------
# Numerical audit
# ------------------------------------------------------------

counts = {}

for stress in STRESSES:
    column = f"{stress}_response"

    counts[stress] = {
        "GAIN": int((master[column] == "GAIN").sum()),
        "LOSS": int((master[column] == "LOSS").sum()),
    }

expected_counts = {
    "HEAT": {"GAIN": 493, "LOSS": 505},
    "COLD": {"GAIN": 138, "LOSS": 248},
    "UVC": {"GAIN": 18, "LOSS": 20},
}

if counts != expected_counts:
    raise RuntimeError(
        f"Unexpected switch counts: {counts}"
    )

membership = pd.DataFrame({
    stress: master[f"{stress}_response"].notna()
    for stress in STRESSES
})

pattern_order = [
    ("HEAT",),
    ("COLD",),
    ("UVC",),
    ("HEAT", "COLD"),
    ("HEAT", "UVC"),
    ("COLD", "UVC"),
    ("HEAT", "COLD", "UVC"),
]

intersection_counts = {}

for pattern in pattern_order:
    mask = np.ones(len(master), dtype=bool)

    for stress in STRESSES:
        if stress in pattern:
            mask &= membership[stress].to_numpy()
        else:
            mask &= ~membership[stress].to_numpy()

    intersection_counts[pattern] = int(mask.sum())

expected_intersections = {
    ("HEAT",): 893,
    ("COLD",): 282,
    ("UVC",): 32,
    ("HEAT", "COLD"): 101,
    ("HEAT", "UVC"): 3,
    ("COLD", "UVC"): 2,
    ("HEAT", "COLD", "UVC"): 1,
}

if intersection_counts != expected_intersections:
    raise RuntimeError(
        "Intersection-count mismatch: "
        f"{intersection_counts}"
    )

locus_rows = (
    recurrent[
        [
            "locus",
            "UVC_response",
            "UVC_delta_A",
            "UVC_delta_B",
            "HEAT_response",
            "COLD_response",
            "gene_ID",
            "gene_symbol",
        ]
    ]
    .drop_duplicates(subset=["locus"])
    .copy()
)

if len(locus_rows) != 6:
    raise RuntimeError(
        f"Expected 6 recurrent loci; found {len(locus_rows)}"
    )

gene_labels = (
    recurrent
    .groupby("locus", sort=False)
    .apply(
        lambda frame:
        "/".join(
            frame["gene_symbol"]
            .replace("", np.nan)
            .dropna()
            .astype(str)
            .tolist()
        )
        or
        "/".join(frame["gene_ID"].astype(str).tolist()),
        include_groups=False
    )
    .to_dict()
)

locus_rows["label"] = locus_rows.apply(
    lambda row:
    f"{row['locus']}  {gene_labels[row['locus']]}",
    axis=1
)

response_value = {
    "LOSS": -1,
    "GAIN": 1,
}

heatmap = []

for _, row in locus_rows.iterrows():
    heatmap.append([
        response_value.get(row["HEAT_response"], 0),
        response_value.get(row["COLD_response"], 0),
        response_value.get(row["UVC_response"], 0),
    ])

heatmap = np.asarray(heatmap)


# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

fig = plt.figure(figsize=(13.2, 9.2))

grid = fig.add_gridspec(
    2,
    2,
    height_ratios=[1.0, 1.15],
    width_ratios=[1.0, 1.2],
    hspace=0.42,
    wspace=0.38
)

ax_a = fig.add_subplot(grid[0, 0])
ax_b = fig.add_subplot(grid[0, 1])
ax_c = fig.add_subplot(grid[1, 0])
ax_d = fig.add_subplot(grid[1, 1])


# Panel A: burden
x = np.arange(len(STRESSES))

gain = np.array([
    counts[stress]["GAIN"]
    for stress in STRESSES
])

loss = np.array([
    counts[stress]["LOSS"]
    for stress in STRESSES
])

ax_a.bar(
    x,
    gain,
    color=COLORS["GAIN"],
    label="CG gain"
)

ax_a.bar(
    x,
    loss,
    bottom=gain,
    color=COLORS["LOSS"],
    label="CG loss"
)

for index, stress in enumerate(STRESSES):
    total = gain[index] + loss[index]

    ax_a.text(
        index,
        total + max(total * 0.025, 8),
        f"{total:,}",
        ha="center",
        va="bottom",
        fontsize=9
    )

ax_a.set_xticks(x)
ax_a.set_xticklabels(["Heat", "Cold", "UVC"])
ax_a.set_ylabel("Stress-specific 100-bp CG switches")
ax_a.set_title("A  Stress-specific CG-switch burden", loc="left")
ax_a.legend(frameon=False)
ax_a.spines[["top", "right"]].set_visible(False)


# Panel B: exact intersections
pattern_labels = [
    "Heat only",
    "Cold only",
    "UVC only",
    "Heat + cold",
    "Heat + UVC",
    "Cold + UVC",
    "All three",
]

values = [
    intersection_counts[pattern]
    for pattern in pattern_order
]

bar_colors = [
    "#D9A441",
    "#56B4E9",
    "#009E73",
    "#9C755F",
    "#CC79A7",
    "#6C8EBF",
    "#332288",
]

bars = ax_b.bar(
    np.arange(len(values)),
    values,
    color=bar_colors
)

ax_b.set_yscale("log")
ax_b.set_ylabel("Exact 100-bp loci (log scale)")
ax_b.set_xticks(np.arange(len(values)))
ax_b.set_xticklabels(
    pattern_labels,
    rotation=35,
    ha="right"
)
ax_b.set_title("B  Exact cross-stress intersections", loc="left")
ax_b.spines[["top", "right"]].set_visible(False)

for bar, value in zip(bars, values):
    ax_b.text(
        bar.get_x() + bar.get_width() / 2,
        value * 1.18,
        str(value),
        ha="center",
        va="bottom",
        fontsize=8
    )


# Panel C: observed versus expected
comparison_labels = [
    "Heat–cold",
    "Heat–UVC",
    "Cold–UVC",
]

observed = enrich["observed_shared_switches"].to_numpy(
    dtype=float
)

expected = enrich["expected_shared_switches"].to_numpy(
    dtype=float
)

position = np.arange(len(comparison_labels))
width = 0.36

ax_c.bar(
    position - width / 2,
    observed,
    width,
    color="#4C78A8",
    label="Observed"
)

ax_c.bar(
    position + width / 2,
    expected,
    width,
    color="#BDBDBD",
    label="Expected"
)

ax_c.set_yscale("log")
ax_c.set_xticks(position)
ax_c.set_xticklabels(comparison_labels)
ax_c.set_ylabel("Shared switches (log scale)")
ax_c.set_title("C  Overlap enrichment", loc="left")
ax_c.legend(frameon=False)
ax_c.spines[["top", "right"]].set_visible(False)

for index, row in enrich.iterrows():
    ax_c.text(
        index,
        observed[index] * 1.25,
        (
            f"{row['fold_enrichment']:.0f}×\n"
            f"P={row['Fisher_upper_P']:.1e}"
        ),
        ha="center",
        va="bottom",
        fontsize=7.5
    )


# Panel D: recurrent loci
cmap = ListedColormap([
    COLORS["LOSS"],
    COLORS["ABSENT"],
    COLORS["GAIN"],
])

image = ax_d.imshow(
    heatmap,
    cmap=cmap,
    vmin=-1,
    vmax=1,
    aspect="auto",
    interpolation="nearest"
)

ax_d.set_xticks(np.arange(3))
ax_d.set_xticklabels(["Heat", "Cold", "UVC"])

ax_d.set_yticks(np.arange(len(locus_rows)))
ax_d.set_yticklabels(locus_rows["label"])

ax_d.set_title(
    "D  Recurrent stress-responsive CG loci",
    loc="left"
)

for row_index in range(heatmap.shape[0]):
    for column_index in range(heatmap.shape[1]):
        value = heatmap[row_index, column_index]

        label = {
            -1: "Loss",
            0: "—",
            1: "Gain",
        }[value]

        color = "white" if value != 0 else "#555555"

        ax_d.text(
            column_index,
            row_index,
            label,
            ha="center",
            va="center",
            color=color,
            fontsize=8,
            fontweight="bold" if value != 0 else "normal"
        )

ax_d.set_xticks(
    np.arange(-0.5, 3, 1),
    minor=True
)
ax_d.set_yticks(
    np.arange(-0.5, len(locus_rows), 1),
    minor=True
)
ax_d.grid(
    which="minor",
    color="white",
    linewidth=1.5
)
ax_d.tick_params(which="minor", bottom=False, left=False)

legend_elements = [
    Line2D(
        [0], [0],
        marker="s",
        color="none",
        markerfacecolor=COLORS["GAIN"],
        markeredgecolor="none",
        markersize=9,
        label="Gain"
    ),
    Line2D(
        [0], [0],
        marker="s",
        color="none",
        markerfacecolor=COLORS["LOSS"],
        markeredgecolor="none",
        markersize=9,
        label="Loss"
    ),
    Line2D(
        [0], [0],
        marker="s",
        color="none",
        markerfacecolor=COLORS["ABSENT"],
        markeredgecolor="none",
        markersize=9,
        label="No switch"
    ),
]

ax_d.legend(
    handles=legend_elements,
    frameon=False,
    ncol=3,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.13)
)


fig.suptitle(
    "Conserved and stress-specific CG methylation responses "
    "after 25 generations",
    fontsize=14,
    y=0.985
)

fig.savefig(
    PNG,
    dpi=400,
    bbox_inches="tight"
)

fig.savefig(
    PDF,
    bbox_inches="tight"
)

plt.close(fig)


# ------------------------------------------------------------
# Source-data audit
# ------------------------------------------------------------

with AUDIT.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow([
        "panel",
        "metric",
        "value",
        "source"
    ])

    for stress in STRESSES:
        for response in ["GAIN", "LOSS"]:
            writer.writerow([
                "A",
                f"{stress}_{response}",
                counts[stress][response],
                MASTER
            ])

    for pattern in pattern_order:
        writer.writerow([
            "B",
            "+".join(pattern),
            intersection_counts[pattern],
            MASTER
        ])

    for _, row in enrich.iterrows():
        writer.writerow([
            "C",
            (
                f"{row['comparison']}_"
                "observed_expected_fold_P"
            ),
            (
                f"{row['observed_shared_switches']};"
                f"{row['expected_shared_switches']};"
                f"{row['fold_enrichment']};"
                f"{row['Fisher_upper_P']}"
            ),
            ENRICH
        ])

    writer.writerow([
        "D",
        "recurrent_loci",
        len(locus_rows),
        RECURRENT
    ])


print("===== NUMERICAL AUDIT =====")

for stress in STRESSES:
    print(
        stress,
        f"GAIN={counts[stress]['GAIN']}",
        f"LOSS={counts[stress]['LOSS']}",
        f"TOTAL={sum(counts[stress].values())}"
    )

print("\nExact intersection counts:")

for pattern in pattern_order:
    print(
        "+".join(pattern),
        intersection_counts[pattern]
    )

print(f"\nRecurrent loci: {len(locus_rows)}")
print(f"WROTE: {PNG}")
print(f"WROTE: {PDF}")
print(f"WROTE: {AUDIT}")
print("CROSS-STRESS CG FIGURE: COMPLETE")
