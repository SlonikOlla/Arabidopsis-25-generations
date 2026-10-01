#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import unquote
import csv

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")

BASE = ROOT / "results/cross_stress_CG_comparison"

INPUT = (
    BASE /
    "CG_recurrent_UVC_heat_cold_loci_annotation.tsv"
)

GFF = (
    ROOT /
    "reference/TAIR10/annotation/"
    "Araport11_GFF3_genes_transposons.20240701.gff"
)

OUTPUT = (
    BASE /
    "CG_recurrent_UVC_heat_cold_gene_descriptions.tsv"
)


def parse_attributes(text):
    attributes = {}

    for item in text.strip().split(";"):
        if "=" not in item:
            continue

        key, value = item.split("=", 1)
        attributes[key] = unquote(value)

    return attributes


recurrent = []

with INPUT.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        gene_ids = [
            gene_id
            for gene_id in row["gene_IDs"].split(",")
            if gene_id
        ]

        for gene_id in gene_ids:
            recurrent.append({
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "locus": row["locus"],
                "UVC_response": row["UVC_response"],
                "UVC_delta_A": row["UVC_delta_A"],
                "UVC_delta_B": row["UVC_delta_B"],
                "HEAT_response": row["HEAT_response"],
                "COLD_response": row["COLD_response"],
                "gene_ID": gene_id,
            })

target_genes = {
    row["gene_ID"]
    for row in recurrent
}

annotations = {}

with GFF.open(encoding="latin-1") as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        feature_type = fields[2]

        if feature_type not in {
            "gene",
            "transposable_element_gene"
        }:
            continue

        attributes = parse_attributes(fields[8])
        gene_id = attributes.get("ID", "")

        if gene_id not in target_genes:
            continue

        annotations[gene_id] = {
            "Araport_feature_type": feature_type,
            "gene_symbol": attributes.get(
                "symbol", ""
            ),
            "full_name": attributes.get(
                "full_name", ""
            ),
            "computational_description": attributes.get(
                "computational_description", ""
            ),
            "curator_summary": attributes.get(
                "curator_summary", ""
            ),
            "note": attributes.get(
                "Note", ""
            ),
            "locus_type": attributes.get(
                "locus_type", ""
            ),
            "strand": fields[6],
            "gene_start": fields[3],
            "gene_end": fields[4],
        }

missing = sorted(target_genes - set(annotations))

if missing:
    raise RuntimeError(
        "Missing Araport annotations: " +
        ", ".join(missing)
    )

fields = [
    "chrom",
    "start",
    "end",
    "locus",
    "UVC_response",
    "UVC_delta_A",
    "UVC_delta_B",
    "HEAT_response",
    "COLD_response",
    "gene_ID",
    "gene_symbol",
    "full_name",
    "computational_description",
    "curator_summary",
    "note",
    "locus_type",
    "Araport_feature_type",
    "strand",
    "gene_start",
    "gene_end",
]

with OUTPUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()

    for row in recurrent:
        output = dict(row)
        output.update(annotations[row["gene_ID"]])
        writer.writerow(output)

print("===== RECURRENT LOCUS–GENE RECORDS =====")
print(f"Loci: {len({row['locus'] for row in recurrent})}")
print(f"Genes: {len(target_genes)}")
print(f"Output records: {len(recurrent)}")

print("\n===== FUNCTIONAL ANNOTATIONS =====")

for row in recurrent:
    annotation = annotations[row["gene_ID"]]

    description = (
        annotation["full_name"] or
        annotation["curator_summary"] or
        annotation["computational_description"] or
        annotation["note"] or
        "No functional description"
    )

    print(
        row["locus"],
        row["UVC_response"],
        f"heat={row['HEAT_response'] or '-'}",
        f"cold={row['COLD_response'] or '-'}",
        row["gene_ID"],
        annotation["gene_symbol"] or "-",
        annotation["locus_type"] or "-",
        description,
        sep="\t"
    )

print(f"\nWROTE: {OUTPUT}")
print("RECURRENT STRESS-GENE DESCRIPTION: COMPLETE")
