#!/usr/bin/env python3

from pathlib import Path
from collections import Counter, defaultdict
from statistics import mean, median, stdev
import csv

from scipy.stats import fisher_exact, kruskal

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
INDIR = ROOT / "results/wgs_variants/PRIMARY23_mutation_analysis"

VARIANTS = INDIR / "PRIMARY23_strict_sample_private_homozygous_variants.tsv"
SAMPLES = INDIR / "PRIMARY23_strict_sample_private_summary.tsv"

GROUP_STATS = INDIR / "PRIMARY23_strict_private_group_statistics.tsv"
SPECTRUM = INDIR / "PRIMARY23_strict_private_SNP_spectrum_by_group.tsv"
TESTS = INDIR / "PRIMARY23_strict_private_statistical_tests.tsv"

group_order = ["F2C", "F25C", "F25H", "F25Cd", "F25UV"]
classes = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

def complement(base):
    return {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }[base]

sample_rows = []

with SAMPLES.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        for column in [
            "total", "SNP", "insertion", "deletion", "complex"
        ]:
            row[column] = int(row[column])

        sample_rows.append(row)

group_values = defaultdict(lambda: defaultdict(list))

for row in sample_rows:
    group = row["group"]

    for column in ["total", "SNP", "insertion", "deletion"]:
        group_values[group][column].append(row[column])

with GROUP_STATS.open("w") as handle:
    handle.write(
        "group\tn\tmetric\tmean\tSD\tmedian\tminimum\tmaximum\n"
    )

    for group in group_order:
        for metric in ["total", "SNP", "insertion", "deletion"]:
            values = group_values[group][metric]
            sd = stdev(values) if len(values) > 1 else 0.0

            handle.write(
                f"{group}\t{len(values)}\t{metric}\t"
                f"{mean(values):.4f}\t{sd:.4f}\t"
                f"{median(values):.4f}\t"
                f"{min(values)}\t{max(values)}\n"
            )

spectra = defaultdict(Counter)

with VARIANTS.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        if row["variant_type"] != "SNP":
            continue

        ref = row["ref"]
        alt = row["alt"]

        if ref in {"A", "G"}:
            ref = complement(ref)
            alt = complement(alt)

        mutation_class = f"{ref}>{alt}"
        spectra[row["group"]][mutation_class] += 1

with SPECTRUM.open("w") as handle:
    handle.write(
        "group\tmutation_class\tcount\tgroup_SNP_total\tpercent\n"
    )

    for group in group_order:
        total = sum(spectra[group].values())

        for mutation_class in classes:
            count = spectra[group][mutation_class]
            percent = 100.0 * count / total if total else 0.0

            handle.write(
                f"{group}\t{mutation_class}\t{count}\t"
                f"{total}\t{percent:.4f}\n"
            )

uvc_ct = spectra["F25UV"]["C>T"]
uvc_other = sum(spectra["F25UV"].values()) - uvc_ct

nonuvc_ct = sum(
    spectra[group]["C>T"]
    for group in group_order
    if group != "F25UV"
)
nonuvc_total = sum(
    sum(spectra[group].values())
    for group in group_order
    if group != "F25UV"
)
nonuvc_other = nonuvc_total - nonuvc_ct

odds_ratio, fisher_p = fisher_exact(
    [[uvc_ct, uvc_other], [nonuvc_ct, nonuvc_other]],
    alternative="two-sided",
)

snp_groups = [
    group_values[group]["SNP"]
    for group in group_order
]

kruskal_result = kruskal(*snp_groups)

with TESTS.open("w") as handle:
    handle.write("test\tstatistic\tP_value\tdetails\n")
    handle.write(
        f"UVC_C_to_T_vs_nonUVC_Fisher\t"
        f"{odds_ratio:.8g}\t{fisher_p:.8g}\t"
        f"UVC={uvc_ct}/{uvc_ct + uvc_other};"
        f"nonUVC={nonuvc_ct}/{nonuvc_total}\n"
    )
    handle.write(
        f"Private_SNP_burden_Kruskal_Wallis\t"
        f"{kruskal_result.statistic:.8g}\t"
        f"{kruskal_result.pvalue:.8g}\t"
        f"groups={','.join(group_order)}\n"
    )

print(f"WROTE: {GROUP_STATS}")
print(f"WROTE: {SPECTRUM}")
print(f"WROTE: {TESTS}")
