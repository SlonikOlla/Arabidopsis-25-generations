#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import unquote
import csv

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
BASE = ROOT / "results/wgs_variants/PRIMARY23_mutation_analysis"

ANN = BASE / "PRIMARY23_UVC_confirmed_24_Araport11_ANN.tsv"
GFF = ROOT / "reference/TAIR10/annotation/Araport11_GFF3_genes_transposons.20240701.gff"
OUT = BASE / "PRIMARY23_UVC_confirmed_24_prioritized_functional_annotation.tsv"

priority = {
    "stop_gained": 1,
    "stop_lost": 2,
    "start_lost": 3,
    "splice_acceptor_variant": 4,
    "splice_donor_variant": 5,
    "missense_variant": 6,
    "splice_region_variant&intron_variant": 7,
    "synonymous_variant": 8,
    "5_prime_UTR_variant": 9,
    "3_prime_UTR_variant": 10,
    "intron_variant": 11,
    "upstream_gene_variant": 12,
    "downstream_gene_variant": 13,
    "intergenic_region": 14,
    "intragenic_variant": 15,
}

gene_metadata = {}

with GFF.open(encoding="latin-1") as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue

        attributes = {}
        for item in fields[8].split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attributes[key] = unquote(value)

        gene_id = attributes.get("ID", "")
        gene_metadata[gene_id] = {
            "symbol": attributes.get("symbol", ""),
            "full_name": attributes.get("full_name", ""),
            "description": attributes.get(
                "computational_description",
                attributes.get("Note", "")
            ),
            "locus_type": attributes.get("locus_type", ""),
        }

rows = []

with ANN.open() as handle:
    reader = csv.reader(handle, delimiter="\t")

    for record in reader:
        chrom, pos, variant_id, ref, alt, lineage, sample, ann_text = record

        annotations = []

        for entry in ann_text.split(","):
            field = entry.split("|") + [""] * 16

            effect = field[1]
            distance_text = field[14]
            distance = (
                int(distance_text)
                if distance_text.isdigit()
                else 10**12
            )

            annotations.append({
                "effect": effect,
                "impact": field[2],
                "gene_name": field[3],
                "gene_id": field[4],
                "feature_type": field[5],
                "transcript": field[6],
                "biotype": field[7],
                "hgvs_c": field[9],
                "hgvs_p": field[10],
                "distance": distance_text,
                "warnings": field[15],
                "rank": priority.get(effect, 100),
                "sort_distance": distance,
            })

        primary = min(
            annotations,
            key=lambda item: (
                item["rank"],
                item["sort_distance"],
                item["gene_id"],
                item["transcript"],
            ),
        )

        metadata = gene_metadata.get(primary["gene_id"], {})

        all_effects = sorted({
            item["effect"] for item in annotations
            if item["effect"] != "intragenic_variant"
        })

        all_genes = sorted({
            item["gene_id"] for item in annotations
            if item["gene_id"] and
               "-Protein" not in item["gene_id"] and
               item["effect"] != "intergenic_region"
        })

        rows.append({
            "chrom": chrom,
            "position": pos,
            "variant_ID": variant_id,
            "ref": ref,
            "alt": alt,
            "lineage": lineage,
            "sample": sample,
            "primary_effect": primary["effect"],
            "impact": primary["impact"],
            "gene_ID": primary["gene_id"],
            "gene_symbol": metadata.get("symbol", ""),
            "gene_full_name": metadata.get("full_name", ""),
            "gene_description": metadata.get("description", ""),
            "locus_type": metadata.get("locus_type", ""),
            "transcript": primary["transcript"],
            "biotype": primary["biotype"],
            "HGVS_c": primary["hgvs_c"],
            "HGVS_p": primary["hgvs_p"],
            "distance_bp": primary["distance"],
            "warnings": primary["warnings"],
            "all_gene_IDs": ",".join(all_genes),
            "all_effects": ",".join(all_effects),
        })

fieldnames = list(rows[0])

with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

if len(rows) != 24:
    raise RuntimeError(f"Expected 24 variants, found {len(rows)}")

print(f"WROTE: {OUT}")
print(f"RECORDS: {len(rows)}")
