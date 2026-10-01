#!/usr/bin/env python3

import csv
import math
from pathlib import Path
from statistics import mean, stdev

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
QC = ROOT / "results/wgs_variants/PRIMARY23_QC_raw"
META = ROOT / "metadata/wgs_primary23_groups.tsv"

MISSING = QC / "PRIMARY23_raw.missingness.imiss"
DEPTH = QC / "PRIMARY23_raw.depth.idepth"
HET = QC / "PRIMARY23_raw.heterozygosity.het"

SAMPLE_OUT = QC / "PRIMARY23_raw_sample_QC_summary.tsv"
GROUP_OUT = QC / "PRIMARY23_raw_group_QC_summary.tsv"

groups = {}

with META.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        groups[row["sample"]] = row["group"]

missingness = {}

with MISSING.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        missingness[row["INDV"]] = {
            "n_data": int(row["N_DATA"]),
            "n_missing": int(row["N_MISS"]),
            "missing_fraction": float(row["F_MISS"]),
        }

depth = {}

with DEPTH.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        depth[row["INDV"]] = {
            "depth_sites": int(row["N_SITES"]),
            "mean_depth": float(row["MEAN_DEPTH"]),
        }

heterozygosity = {}

with HET.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        sample = row["INDV"]
        observed_hom = int(row["O(HOM)"])
        expected_hom = float(row["E(HOM)"])
        n_sites = int(row["N_SITES"])
        f_value = float(row["F"])

        observed_het = n_sites - observed_hom
        observed_het_fraction = (
            observed_het / n_sites if n_sites else math.nan
        )

        heterozygosity[sample] = {
            "het_sites": n_sites,
            "observed_hom": observed_hom,
            "expected_hom": expected_hom,
            "observed_het": observed_het,
            "observed_het_fraction": observed_het_fraction,
            "f_value": f_value,
        }

sample_sets = [
    set(groups),
    set(missingness),
    set(depth),
    set(heterozygosity),
]

if not all(sample_set == sample_sets[0] for sample_set in sample_sets[1:]):
    raise RuntimeError("Sample sets differ among metadata and QC files")

rows = []

for sample in sorted(groups):
    row = {
        "sample": sample,
        "group": groups[sample],
        **missingness[sample],
        **depth[sample],
        **heterozygosity[sample],
    }
    rows.append(row)

sample_fields = [
    "sample",
    "group",
    "n_data",
    "n_missing",
    "missing_fraction",
    "depth_sites",
    "mean_depth",
    "het_sites",
    "observed_hom",
    "expected_hom",
    "observed_het",
    "observed_het_fraction",
    "f_value",
]

with SAMPLE_OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=sample_fields,
    )
    writer.writeheader()

    for row in rows:
        formatted = dict(row)
        formatted["missing_fraction"] = (
            f'{row["missing_fraction"]:.6f}'
        )
        formatted["mean_depth"] = f'{row["mean_depth"]:.3f}'
        formatted["expected_hom"] = f'{row["expected_hom"]:.1f}'
        formatted["observed_het_fraction"] = (
            f'{row["observed_het_fraction"]:.6f}'
        )
        formatted["f_value"] = f'{row["f_value"]:.5f}'
        writer.writerow(formatted)

group_order = ["F2C", "F25C", "F25H", "F25Cd", "F25UV"]
group_fields = [
    "group",
    "n",
    "missing_mean",
    "missing_sd",
    "depth_mean",
    "depth_sd",
    "observed_het_fraction_mean",
    "observed_het_fraction_sd",
]

with GROUP_OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=group_fields,
    )
    writer.writeheader()

    for group in group_order:
        subset = [row for row in rows if row["group"] == group]

        missing_values = [
            row["missing_fraction"] for row in subset
        ]
        depth_values = [
            row["mean_depth"] for row in subset
        ]
        het_values = [
            row["observed_het_fraction"] for row in subset
        ]

        writer.writerow({
            "group": group,
            "n": len(subset),
            "missing_mean": f"{mean(missing_values):.6f}",
            "missing_sd": f"{stdev(missing_values):.6f}",
            "depth_mean": f"{mean(depth_values):.3f}",
            "depth_sd": f"{stdev(depth_values):.3f}",
            "observed_het_fraction_mean": (
                f"{mean(het_values):.6f}"
            ),
            "observed_het_fraction_sd": (
                f"{stdev(het_values):.6f}"
            ),
        })

print(f"WROTE: {SAMPLE_OUT}")
print(f"WROTE: {GROUP_OUT}")
print(f"SAMPLES: {len(rows)}")
print("PRIMARY23 RAW QC SUMMARY: COMPLETE")
