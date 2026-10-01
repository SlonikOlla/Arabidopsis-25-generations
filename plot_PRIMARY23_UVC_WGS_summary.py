#!/usr/bin/env python3

from pathlib import Path
import itertools
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
BASE = ROOT / "results/wgs_variants/PRIMARY23_mutation_analysis"
FIGDIR = BASE / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

SAMPLE_FILE = BASE / "PRIMARY23_UVC_confirmed_SNP_sample_summary.tsv"
SPECTRUM_FILE = BASE / "PRIMARY23_UVC_confirmed_SNP_spectrum.tsv"
FUNCTION_FILE = BASE / "PRIMARY23_UVC_confirmed_24_final_functional_annotation.tsv"
CONTEXT_FILE = BASE / "PRIMARY23_UVC_confirmed_SNP_context.tsv"

PNG = FIGDIR / "PRIMARY23_UVC_WGS_summary.png"
PDF = FIGDIR / "PRIMARY23_UVC_WGS_summary.pdf"
SOURCE = BASE / "PRIMARY23_UVC_WGS_figure_source_audit.tsv"

for path in [SAMPLE_FILE, SPECTRUM_FILE, FUNCTION_FILE, CONTEXT_FILE]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty input: {path}")

samples = pd.read_csv(SAMPLE_FILE, sep="\t")
spectrum = pd.read_csv(SPECTRUM_FILE, sep="\t")
functional = pd.read_csv(FUNCTION_FILE, sep="\t")
context = pd.read_csv(CONTEXT_FILE, sep="\t")

required_sample = {"sample", "lineage", "confirmed_SNPs"}
required_spectrum = {"mutation_class", "count", "percent"}
required_functional = {"variant_ID", "corrected_primary_effect"}
required_context = {
    "chrom", "position", "sample",
    "normalized_trinucleotide", "mutation_class"
}

for observed, required, name in [
    (set(samples.columns), required_sample, "sample table"),
    (set(spectrum.columns), required_spectrum, "spectrum table"),
    (set(functional.columns), required_functional, "functional table"),
    (set(context.columns), required_context, "context table"),
]:
    missing = required - observed
    if missing:
        raise RuntimeError(f"{name} missing columns: {sorted(missing)}")

sample_order = [
    "F25UV_7_WGS",
    "F25UV_8_WGS",
    "F25UV_10_WGS",
    "F25UV_11_WGS",
    "F25UV_13_WGS",
]

samples = samples.set_index("sample").loc[sample_order].reset_index()

spectrum_order = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
spectrum = (
    spectrum.set_index("mutation_class")
    .reindex(spectrum_order)
    .reset_index()
)

effect_order = [
    "transposable_element_exonic_variant",
    "upstream_gene_variant",
    "missense_variant",
    "intron_variant",
    "downstream_gene_variant",
    "synonymous_variant",
    "splice_region_variant&intron_variant",
    "intergenic_region",
    "5_prime_UTR_variant",
]

effect_counts = (
    functional["corrected_primary_effect"]
    .value_counts()
    .reindex(effect_order, fill_value=0)
)

effect_labels = {
    "transposable_element_exonic_variant": "TE-gene exonic",
    "upstream_gene_variant": "Upstream",
    "missense_variant": "Missense",
    "intron_variant": "Intronic",
    "downstream_gene_variant": "Downstream",
    "synonymous_variant": "Synonymous",
    "splice_region_variant&intron_variant": "Splice-region/intronic",
    "intergenic_region": "Intergenic",
    "5_prime_UTR_variant": "5′ UTR",
}

ct_context = context.loc[
    context["mutation_class"] == "C>T",
    "normalized_trinucleotide"
]

canonical_contexts = [
    left + "C" + right
    for left, right in itertools.product("ACGT", repeat=2)
]

context_counts = (
    ct_context.value_counts()
    .reindex(canonical_contexts, fill_value=0)
)

# Numerical validation
if len(samples) != 5:
    raise RuntimeError(f"Expected 5 UVC samples, found {len(samples)}")

if samples["confirmed_SNPs"].sum() != 24:
    raise RuntimeError(
        f"Sample counts sum to {samples['confirmed_SNPs'].sum()}, expected 24"
    )

if spectrum["count"].sum() != 24:
    raise RuntimeError(
        f"Spectrum counts sum to {spectrum['count'].sum()}, expected 24"
    )

if len(functional) != 24 or effect_counts.sum() != 24:
    raise RuntimeError("Functional-annotation count does not equal 24")

if len(context) != 24:
    raise RuntimeError(f"Expected 24 context records, found {len(context)}")

if context_counts.sum() != 18:
    raise RuntimeError(
        f"C>T context counts sum to {context_counts.sum()}, expected 18"
    )

# Frozen figure-source audit
audit_rows = []

for _, row in samples.iterrows():
    audit_rows.append({
        "panel": "A",
        "category": row["sample"],
        "count": int(row["confirmed_SNPs"]),
    })

for _, row in spectrum.iterrows():
    audit_rows.append({
        "panel": "B",
        "category": row["mutation_class"],
        "count": int(row["count"]),
    })

for effect, count in effect_counts.items():
    audit_rows.append({
        "panel": "C",
        "category": effect,
        "count": int(count),
    })

for tri, count in context_counts.items():
    audit_rows.append({
        "panel": "D",
        "category": tri,
        "count": int(count),
    })

pd.DataFrame(audit_rows).to_csv(SOURCE, sep="\t", index=False)

# Plot
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

fig, axes = plt.subplots(
    2, 2,
    figsize=(12.0, 9.5),
    constrained_layout=True
)

ax_a, ax_b, ax_c, ax_d = axes.flat

# Panel A
lineage_colors = {
    "UVC_A": "#377EB8",
    "UVC_B": "#E41A1C",
}

sample_labels = [
    name.replace("F25UV_", "").replace("_WGS", "")
    for name in samples["sample"]
]

bars = ax_a.bar(
    sample_labels,
    samples["confirmed_SNPs"],
    color=[lineage_colors[x] for x in samples["lineage"]],
    edgecolor="black",
    linewidth=0.7
)

for bar, value in zip(bars, samples["confirmed_SNPs"]):
    ax_a.text(
        bar.get_x() + bar.get_width()/2,
        value + 0.25,
        str(int(value)),
        ha="center",
        va="bottom",
        fontsize=9
    )

ax_a.set_ylabel("Confirmed SNPs")
ax_a.set_xlabel("F25UV sample")
ax_a.set_title("A  Independently confirmed SNP burden", loc="left", weight="bold")
ax_a.set_ylim(0, max(samples["confirmed_SNPs"]) + 2)
ax_a.spines[["top", "right"]].set_visible(False)

# Panel B
spectrum_colors = [
    "#4C78A8", "#72B7B2", "#E45756",
    "#F2CF5B", "#54A24B", "#B279A2"
]

bars = ax_b.bar(
    spectrum["mutation_class"],
    spectrum["count"],
    color=spectrum_colors,
    edgecolor="black",
    linewidth=0.7
)

for bar, count, percent in zip(
    bars, spectrum["count"], spectrum["percent"]
):
    ax_b.text(
        bar.get_x() + bar.get_width()/2,
        count + 0.35,
        f"{int(count)}\n({percent:.1f}%)",
        ha="center",
        va="bottom",
        fontsize=8
    )

ax_b.set_ylabel("Confirmed SNPs")
ax_b.set_xlabel("Mutation class")
ax_b.set_title("B  Six-class mutation spectrum", loc="left", weight="bold")
ax_b.set_ylim(0, max(spectrum["count"]) + 4)
ax_b.spines[["top", "right"]].set_visible(False)

# Panel C
display_effects = [effect_labels[x] for x in effect_order]
y_positions = range(len(effect_order))

bars = ax_c.barh(
    list(y_positions),
    effect_counts.values,
    color="#7A5195",
    edgecolor="black",
    linewidth=0.6
)

ax_c.set_yticks(list(y_positions))
ax_c.set_yticklabels(display_effects)
ax_c.invert_yaxis()

for bar, value in zip(bars, effect_counts.values):
    ax_c.text(
        value + 0.10,
        bar.get_y() + bar.get_height()/2,
        str(int(value)),
        va="center",
        fontsize=9
    )

ax_c.set_xlabel("Confirmed SNPs")
ax_c.set_title("C  Corrected functional consequences", loc="left", weight="bold")
ax_c.set_xlim(0, max(effect_counts.values) + 1.3)
ax_c.spines[["top", "right"]].set_visible(False)

# Panel D
context_colors = [
    "#F4A582" if count > 0 else "#E0E0E0"
    for count in context_counts.values
]

bars = ax_d.bar(
    canonical_contexts,
    context_counts.values,
    color=context_colors,
    edgecolor="black",
    linewidth=0.5
)

for bar, value in zip(bars, context_counts.values):
    if value > 0:
        ax_d.text(
            bar.get_x() + bar.get_width()/2,
            value + 0.08,
            str(int(value)),
            ha="center",
            va="bottom",
            fontsize=8
        )

ax_d.set_ylabel("C>T substitutions")
ax_d.set_xlabel("Pyrimidine-oriented trinucleotide")
ax_d.set_title("D  C>T sequence context", loc="left", weight="bold")
ax_d.tick_params(axis="x", rotation=45)
ax_d.set_ylim(0, max(context_counts.values) + 1)
ax_d.spines[["top", "right"]].set_visible(False)

fig.savefig(PNG, dpi=600, bbox_inches="tight")
fig.savefig(PDF, bbox_inches="tight")
plt.close(fig)

print("===== NUMERICAL AUDIT =====")
print(f"Sample records: {len(samples)}")
print(f"Confirmed SNPs across samples: {samples['confirmed_SNPs'].sum()}")
print(f"Spectrum total: {spectrum['count'].sum()}")
print(f"Functional records: {len(functional)}")
print(f"C>T records: {context_counts.sum()}")
print(f"Nonzero C>T contexts: {(context_counts > 0).sum()}")

print("\n===== SAMPLE COUNTS =====")
for _, row in samples.iterrows():
    print(row["sample"], row["lineage"], int(row["confirmed_SNPs"]))

print("\n===== FUNCTIONAL COUNTS =====")
for effect, count in effect_counts.items():
    print(effect, int(count))

print("\n===== C>T CONTEXT COUNTS =====")
for tri, count in context_counts.items():
    print(tri, int(count))

print("\n===== OUTPUTS =====")
print(PNG)
print(PDF)
print(SOURCE)
print("PRIMARY23 UVC WGS SUMMARY FIGURE: COMPLETE")
