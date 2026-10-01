#!/usr/bin/env python3

import csv
from pathlib import Path

root = Path("/scratch/vasilisa/arabidopsis_25gen")
indir = root / "results/wgs_variants/PRIMARY23_mutation_analysis"
src = indir / "PRIMARY23_UVC_matched_genotypes.tsv"
out = indir / "PRIMARY23_UVC_matched_candidate_variants.tsv"

uvc_a = ["F25UV_7_WGS", "F25UV_8_WGS"]
uvc_b = ["F25UV_10_WGS", "F25UV_11_WGS", "F25UV_13_WGS"]

backgrounds = {
    "UVC_A": ["SRR17867643"],
    "UVC_B": ["SRR17867641", "SRR17867645"],
}

def homozygous_allele(row, sample):
    gt = row[f"{sample}_GT"]
    if gt not in {"0/0", "0|0", "1/1", "1|1"}:
        return None

    try:
        dp = int(row[f"{sample}_DP"])
        gq = int(row[f"{sample}_GQ"])
        ad = [int(x) for x in row[f"{sample}_AD"].split(",")[:2]]
    except (ValueError, IndexError):
        return None

    if dp < 8 or gq < 20 or len(ad) < 2:
        return None

    total = ad[0] + ad[1]
    if total == 0:
        return None

    allele = 0 if gt[0] == "0" else 1
    supporting = ad[allele] / total

    if supporting < 0.90:
        return None

    return allele

def emit(writer, row, class_name, lineage, targets, comparators):
    target_calls = [homozygous_allele(row, s) for s in targets]
    comparator_calls = [homozygous_allele(row, s) for s in comparators]

    if None in target_calls or None in comparator_calls:
        return

    if len(set(target_calls)) != 1 or len(set(comparator_calls)) != 1:
        return

    target_allele = target_calls[0]
    comparator_allele = comparator_calls[0]

    if target_allele == comparator_allele:
        return

    direction = (
        "REF_to_ALT" if comparator_allele == 0
        else "ALT_to_REF"
    )

    writer.writerow({
        "chrom": row["chrom"],
        "position": row["position"],
        "ref": row["ref"],
        "alt": row["alt"],
        "type": row["type"],
        "class": class_name,
        "lineage": lineage,
        "target_samples": ",".join(targets),
        "comparator_samples": ",".join(comparators),
        "comparator_allele": comparator_allele,
        "uvc_allele": target_allele,
        "direction": direction,
    })

fields = [
    "chrom", "position", "ref", "alt", "type",
    "class", "lineage", "target_samples", "comparator_samples",
    "comparator_allele", "uvc_allele", "direction",
]

with src.open() as handle, out.open("w", newline="") as output:
    reader = csv.DictReader(handle, delimiter="\t")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for row in reader:
        emit(
            writer, row, "lineage_shared", "UVC_A",
            uvc_a, backgrounds["UVC_A"]
        )
        emit(
            writer, row, "lineage_shared", "UVC_B",
            uvc_b, backgrounds["UVC_B"]
        )

        for sample in uvc_a:
            comparators = [
                s for s in uvc_a if s != sample
            ] + backgrounds["UVC_A"]
            emit(
                writer, row, "sample_specific", "UVC_A",
                [sample], comparators
            )

        for sample in uvc_b:
            comparators = [
                s for s in uvc_b if s != sample
            ] + backgrounds["UVC_B"]
            emit(
                writer, row, "sample_specific", "UVC_B",
                [sample], comparators
            )

print(f"WROTE: {out}")
