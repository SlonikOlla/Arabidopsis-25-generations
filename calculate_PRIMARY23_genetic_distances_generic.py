#!/usr/bin/env python3

import csv
import math
import sys
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: calculate_PRIMARY23_genetic_distances.py "
        "PLINK_RAW GROUPS_TSV OUTPUT_DIR OUTPUT_PREFIX"
    )

raw_file = Path(sys.argv[1])
group_file = Path(sys.argv[2])
output_dir = Path(sys.argv[3])
prefix = sys.argv[4]

output_dir.mkdir(parents=True, exist_ok=True)

metadata_columns = {
    "#FID", "FID", "IID", "PAT", "MAT",
    "SEX", "PHENOTYPE", "PHENO1"
}

groups = {}

with group_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        groups[row["sample"]] = row["group"]

with raw_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    header = reader.fieldnames

    if header is None:
        raise RuntimeError("PLINK raw file has no header")

    if "IID" not in header:
        raise RuntimeError("IID column absent from PLINK raw file")

    genotype_columns = [
        column for column in header
        if column not in metadata_columns
    ]

    samples = []
    genotypes = {}

    for row in reader:
        sample = row["IID"]

        values = []

        for column in genotype_columns:
            value = row[column]

            if value in {"", "NA", "nan", "."}:
                values.append(None)
            else:
                values.append(float(value))

        samples.append(sample)
        genotypes[sample] = values

if len(samples) != 23:
    raise RuntimeError(
        f"Expected 23 samples but found {len(samples)}"
    )

if len(genotype_columns) < 1:
    raise RuntimeError("No genotype columns found")

missing_groups = [
    sample for sample in samples
    if sample not in groups
]

if missing_groups:
    raise RuntimeError(
        "Samples absent from group table: " +
        ", ".join(missing_groups)
    )

pairwise_rows = []
distance_lookup = {
    sample: {} for sample in samples
}

for index_a, sample_a in enumerate(samples):
    genotype_a = genotypes[sample_a]

    for index_b in range(index_a + 1, len(samples)):
        sample_b = samples[index_b]
        genotype_b = genotypes[sample_b]

        paired = [
            (a, b)
            for a, b in zip(genotype_a, genotype_b)
            if a is not None and b is not None
        ]

        n_sites = len(paired)

        if n_sites == 0:
            raise RuntimeError(
                f"No comparable genotypes for {sample_a} and {sample_b}"
            )

        absolute_differences = [
            abs(a - b) for a, b in paired
        ]

        ibs_distance = sum(
            difference / 2.0
            for difference in absolute_differences
        ) / n_sites

        exact_match_fraction = sum(
            1 for difference in absolute_differences
            if difference == 0
        ) / n_sites

        genotype_discordance = sum(
            1 for difference in absolute_differences
            if difference != 0
        ) / n_sites

        mean_a = sum(a for a, _ in paired) / n_sites
        mean_b = sum(b for _, b in paired) / n_sites

        numerator = sum(
            (a - mean_a) * (b - mean_b)
            for a, b in paired
        )

        denominator_a = sum(
            (a - mean_a) ** 2
            for a, _ in paired
        )

        denominator_b = sum(
            (b - mean_b) ** 2
            for _, b in paired
        )

        denominator = math.sqrt(
            denominator_a * denominator_b
        )

        dosage_correlation = (
            numerator / denominator
            if denominator > 0
            else float("nan")
        )

        distance_lookup[sample_a][sample_b] = ibs_distance
        distance_lookup[sample_b][sample_a] = ibs_distance

        pairwise_rows.append({
            "sample_1": sample_a,
            "group_1": groups[sample_a],
            "sample_2": sample_b,
            "group_2": groups[sample_b],
            "joint_sites": n_sites,
            "ibs_distance": ibs_distance,
            "exact_match_fraction": exact_match_fraction,
            "genotype_discordance": genotype_discordance,
            "dosage_correlation": dosage_correlation,
        })

long_output = output_dir / f"{prefix}_pairwise_distances.tsv"

with long_output.open("w", newline="") as handle:
    fieldnames = [
        "sample_1",
        "group_1",
        "sample_2",
        "group_2",
        "joint_sites",
        "ibs_distance",
        "exact_match_fraction",
        "genotype_discordance",
        "dosage_correlation",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()

    for row in sorted(
        pairwise_rows,
        key=lambda value: value["ibs_distance"]
    ):
        output_row = dict(row)

        for field in [
            "ibs_distance",
            "exact_match_fraction",
            "genotype_discordance",
            "dosage_correlation",
        ]:
            output_row[field] = f"{output_row[field]:.8f}"

        writer.writerow(output_row)

matrix_output = output_dir / f"{prefix}_IBS_distance_matrix.tsv"

with matrix_output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample"] + samples)

    for sample_a in samples:
        row = [sample_a]

        for sample_b in samples:
            if sample_a == sample_b:
                row.append("0.00000000")
            else:
                row.append(
                    f"{distance_lookup[sample_a][sample_b]:.8f}"
                )

        writer.writerow(row)

nearest_output = output_dir / f"{prefix}_nearest_neighbors.tsv"

with nearest_output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow([
        "sample",
        "group",
        "nearest_sample",
        "nearest_group",
        "ibs_distance",
        "same_group",
        "second_nearest_sample",
        "second_nearest_group",
        "second_nearest_distance",
    ])

    for sample in samples:
        neighbors = sorted(
            distance_lookup[sample].items(),
            key=lambda value: value[1]
        )

        nearest_sample, nearest_distance = neighbors[0]
        second_sample, second_distance = neighbors[1]

        writer.writerow([
            sample,
            groups[sample],
            nearest_sample,
            groups[nearest_sample],
            f"{nearest_distance:.8f}",
            "YES" if groups[sample] == groups[nearest_sample] else "NO",
            second_sample,
            groups[second_sample],
            f"{second_distance:.8f}",
        ])

uvc_output = output_dir / f"{prefix}_UVC_ranked_neighbors.tsv"

with uvc_output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow([
        "UVC_sample",
        "rank",
        "neighbor",
        "neighbor_group",
        "ibs_distance",
    ])

    for sample in samples:
        if groups[sample] != "F25UV":
            continue

        neighbors = sorted(
            distance_lookup[sample].items(),
            key=lambda value: value[1]
        )

        for rank, (neighbor, distance) in enumerate(
            neighbors,
            start=1
        ):
            writer.writerow([
                sample,
                rank,
                neighbor,
                groups[neighbor],
                f"{distance:.8f}",
            ])

print(f"SAMPLES: {len(samples)}")
print(f"VARIANTS: {len(genotype_columns)}")
print(f"PAIRWISE COMPARISONS: {len(pairwise_rows)}")
print(f"WROTE: {long_output}")
print(f"WROTE: {matrix_output}")
print(f"WROTE: {nearest_output}")
print(f"WROTE: {uvc_output}")
