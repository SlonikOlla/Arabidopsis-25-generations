#!/usr/bin/env python3

from pathlib import Path
from math import comb
import csv

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
FASTA = ROOT / "reference/TAIR10/TAIR10_chr_all.fa"
CONTEXT = (
    ROOT /
    "results/wgs_variants/PRIMARY23_mutation_analysis/"
    "PRIMARY23_UVC_confirmed_SNP_context.tsv"
)
OUTPUT = (
    ROOT /
    "results/wgs_variants/PRIMARY23_mutation_analysis/"
    "PRIMARY23_UVC_dipyrimidine_enrichment.tsv"
)

sequences = {}
name = None
parts = []

with FASTA.open() as handle:
    for line in handle:
        line = line.strip()

        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(parts).upper()

            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)

    if name is not None:
        sequences[name] = "".join(parts).upper()

genomic_c_sites = 0
genomic_dipyrimidine_sites = 0
genomic_cpg_sites = 0

for chromosome in ["1", "2", "3", "4", "5"]:
    sequence = sequences[chromosome]

    for index in range(1, len(sequence) - 1):
        base = sequence[index]

        if base == "C":
            left = sequence[index - 1]
            right = sequence[index + 1]
        elif base == "G":
            left = sequence[index + 1].translate(
                str.maketrans("ACGT", "TGCA")
            )
            right = sequence[index - 1].translate(
                str.maketrans("ACGT", "TGCA")
            )
        else:
            continue

        if left not in "ACGT" or right not in "ACGT":
            continue

        genomic_c_sites += 1

        if left in "CT" or right in "CT":
            genomic_dipyrimidine_sites += 1

        if right == "G":
            genomic_cpg_sites += 1

observed_total = 0
observed_dipyrimidine = 0
observed_cpg = 0

with CONTEXT.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        if row["mutation_class"] != "C>T":
            continue

        observed_total += 1

        if row["dipyrimidine"] == "YES":
            observed_dipyrimidine += 1

        if row["CpG"] == "YES":
            observed_cpg += 1

background_probability = (
    genomic_dipyrimidine_sites / genomic_c_sites
)

binomial_p = sum(
    comb(observed_total, k) *
    background_probability ** k *
    (1.0 - background_probability) ** (observed_total - k)
    for k in range(observed_dipyrimidine, observed_total + 1)
)

observed_fraction = observed_dipyrimidine / observed_total
fold_enrichment = observed_fraction / background_probability

with OUTPUT.open("w") as handle:
    handle.write(
        "metric\tvalue\n"
        f"genomic_C_or_G_sites\t{genomic_c_sites}\n"
        f"genomic_dipyrimidine_sites\t{genomic_dipyrimidine_sites}\n"
        f"genomic_dipyrimidine_fraction\t{background_probability:.8f}\n"
        f"genomic_CpG_sites\t{genomic_cpg_sites}\n"
        f"observed_C_to_T\t{observed_total}\n"
        f"observed_C_to_T_at_dipyrimidines\t{observed_dipyrimidine}\n"
        f"observed_dipyrimidine_fraction\t{observed_fraction:.8f}\n"
        f"dipyrimidine_fold_enrichment\t{fold_enrichment:.6f}\n"
        f"exact_binomial_upper_tail_P\t{binomial_p:.8g}\n"
        f"observed_C_to_T_at_CpG\t{observed_cpg}\n"
    )

print(f"GENOMIC C/G SITES: {genomic_c_sites}")
print(
    "GENOMIC DIPYRIMIDINE FRACTION: "
    f"{background_probability:.4f}"
)
print(
    "OBSERVED C>T AT DIPYRIMIDINES: "
    f"{observed_dipyrimidine}/{observed_total} "
    f"({observed_fraction:.4f})"
)
print(f"FOLD ENRICHMENT: {fold_enrichment:.3f}")
print(f"EXACT BINOMIAL P: {binomial_p:.6g}")
print(f"WROTE: {OUTPUT}")
