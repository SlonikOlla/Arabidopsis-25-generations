#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import re

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
BASE = ROOT / "results/wgs_variants/PRIMARY23_mutation_analysis"

TABLE = BASE / "PRIMARY23_UVC_confirmed_24_prioritized_functional_annotation.tsv"
GFF = ROOT / "reference/TAIR10/annotation/Araport11_GFF3_genes_transposons.20240701.gff"
FASTA = ROOT / "reference/TAIR10/TAIR10_chr_all.fa"
OUT = BASE / "PRIMARY23_UVC_coding_consequence_independent_validation.tsv"

complement = str.maketrans("ACGTN", "TGCAN")

bases = "TCAG"
amino_acids = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)

codon_table = {}
index = 0
for first in bases:
    for second in bases:
        for third in bases:
            codon_table[first + second + third] = amino_acids[index]
            index += 1

genome = {}
name = None
sequence = []

with FASTA.open() as handle:
    for line in handle:
        line = line.rstrip()

        if line.startswith(">"):
            if name is not None:
                genome[name] = "".join(sequence).upper()

            name = line[1:].split()[0]
            sequence = []
        else:
            sequence.append(line)

    if name is not None:
        genome[name] = "".join(sequence).upper()

cds_by_transcript = defaultdict(list)

with GFF.open(encoding="latin-1") as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        field = line.rstrip("\n").split("\t")
        if len(field) < 9 or field[2] != "CDS":
            continue

        chromosome = field[0]
        if chromosome.startswith("Chr"):
            chromosome = chromosome[3:]

        if chromosome == "M":
            chromosome = "Mt"
        elif chromosome == "C":
            chromosome = "Pt"

        attributes = {}
        for item in field[8].split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attributes[key] = value

        for parent in attributes.get("Parent", "").split(","):
            if parent:
                cds_by_transcript[parent].append({
                    "chrom": chromosome,
                    "start": int(field[3]),
                    "end": int(field[4]),
                    "strand": field[6],
                    "phase": field[7],
                })

selected = []

with TABLE.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        if row["primary_effect"] in {
            "missense_variant",
            "synonymous_variant",
        }:
            selected.append(row)

results = []

for row in selected:
    transcript = row["transcript"]
    segments = cds_by_transcript.get(transcript, [])

    if not segments:
        print("SKIP_NO_CDS:", row["variant_ID"], transcript)
        continue

    chromosome = row["chrom"]
    position = int(row["position"])
    ref = row["ref"].upper()
    alt = row["alt"].upper()
    strand = segments[0]["strand"]

    if any(segment["chrom"] != chromosome for segment in segments):
        raise RuntimeError(f"Chromosome mismatch for {transcript}")

    ordered = sorted(
        segments,
        key=lambda segment: segment["start"],
        reverse=(strand == "-"),
    )

    cds_sequence = []
    genomic_positions = []

    for segment in ordered:
        start = segment["start"]
        end = segment["end"]

        if strand == "+":
            positions = list(range(start, end + 1))
            segment_sequence = genome[chromosome][start - 1:end]
        else:
            positions = list(range(end, start - 1, -1))
            segment_sequence = genome[chromosome][start - 1:end]
            segment_sequence = segment_sequence.translate(complement)[::-1]

        cds_sequence.append(segment_sequence)
        genomic_positions.extend(positions)

    cds_sequence = "".join(cds_sequence)

    if position not in genomic_positions:
        raise RuntimeError(
            f"Variant {chromosome}:{position} absent from CDS of {transcript}"
        )

    cds_index = genomic_positions.index(position)
    coding_position = cds_index + 1

    oriented_ref = ref if strand == "+" else ref.translate(complement)
    oriented_alt = alt if strand == "+" else alt.translate(complement)

    observed_ref = cds_sequence[cds_index]

    codon_start = (cds_index // 3) * 3
    ref_codon = cds_sequence[codon_start:codon_start + 3]

    alt_codon_list = list(ref_codon)
    alt_codon_list[cds_index % 3] = oriented_alt
    alt_codon = "".join(alt_codon_list)

    ref_aa = codon_table.get(ref_codon, "X")
    alt_aa = codon_table.get(alt_codon, "X")

    translated = [
        codon_table.get(cds_sequence[i:i + 3], "X")
        for i in range(0, len(cds_sequence) - 2, 3)
    ]

    internal_stops = sum(
        amino_acid == "*"
        for amino_acid in translated[:-1]
    )

    hgvs_match = re.search(r"c\.(\d+)", row["HGVS_c"])
    snpeff_coding_position = (
        int(hgvs_match.group(1))
        if hgvs_match
        else None
    )

    expected_effect = (
        "synonymous_variant"
        if ref_aa == alt_aa
        else "missense_variant"
    )

    status = "PASS"

    if observed_ref != oriented_ref:
        status = "FAIL_REFERENCE"

    if expected_effect != row["primary_effect"]:
        status = "FAIL_EFFECT"

    if (
        snpeff_coding_position is not None and
        snpeff_coding_position != coding_position
    ):
        status = "FAIL_CDS_POSITION"

    results.append({
        "variant_ID": row["variant_ID"],
        "chrom": chromosome,
        "position": position,
        "sample": row["sample"],
        "gene_ID": row["gene_ID"],
        "transcript": transcript,
        "strand": strand,
        "ref": ref,
        "alt": alt,
        "coding_position": coding_position,
        "snpEff_HGVS_c": row["HGVS_c"],
        "ref_codon": ref_codon,
        "alt_codon": alt_codon,
        "ref_AA": ref_aa,
        "alt_AA": alt_aa,
        "independent_effect": expected_effect,
        "snpEff_effect": row["primary_effect"],
        "internal_stop_codons": internal_stops,
        "snpEff_warning": row["warnings"],
        "validation_status": status,
    })

with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(results[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(results)

failed = sum(row["validation_status"] != "PASS" for row in results)

print(f"CODING VARIANTS TESTED: {len(results)}")
print(f"VALIDATION FAILURES: {failed}")
print(f"WROTE: {OUT}")
