#!/usr/bin/env python3

from pathlib import Path
import csv

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")

BASE = ROOT / "results/cross_stress_CG_comparison"
UVC_TABLE = BASE / "CG_UVC_switches_heat_cold_comparison.tsv"

HIGH = (
    ROOT /
    "results/wgbs_PRIMARY22_500bp/"
    "UVC_background_adjusted/"
    "CG_cross_method_high_confidence_DMRs.tsv"
)

FEATURES = (
    ROOT /
    "results/wgbs_PRIMARY22_500bp/"
    "UVC_background_adjusted/"
    "Araport11_gene_and_TE_gene_features.bed"
)

HEAT = (
    ROOT /
    "results/trajectories_100bp_PRIMARY_CLEAN/"
    "CG_HEAT_three_stage_trajectories.tsv"
)

COLD = (
    ROOT /
    "results/trajectories_100bp/"
    "CG_COLD_three_stage_trajectories.tsv"
)

RECURRENT_OUT = (
    BASE /
    "CG_recurrent_UVC_heat_cold_loci_annotation.tsv"
)

HIGH_OUT = (
    BASE /
    "CG_high_confidence_UVC_DMR_historical_switch_audit.tsv"
)


def load_features(path):
    features = []

    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row in reader:
            if not row or row[0].startswith("#"):
                continue

            if len(row) < 9:
                raise RuntimeError(
                    f"Malformed feature row: {row}"
                )

            features.append({
                "chrom": row[0],
                "start": int(row[1]),
                "end": int(row[2]),
                "feature_ID": row[3],
                "feature_type": row[4],
                "strand": row[5],
                "gene_ID": row[6],
                "gene_symbol": row[7],
                "locus_type": row[8],
            })

    return features


def feature_distance(start, end, feature):
    if start < feature["end"] and end > feature["start"]:
        return 0

    if end <= feature["start"]:
        return feature["start"] - end

    return start - feature["end"]


def annotate_interval(chrom, start, end, features):
    chromosome_features = [
        feature for feature in features
        if feature["chrom"] == chrom
    ]

    overlaps = [
        feature for feature in chromosome_features
        if feature_distance(start, end, feature) == 0
    ]

    if overlaps:
        selected = overlaps
        relationship = "DIRECT_OVERLAP"
        distance = 0
    else:
        minimum = min(
            feature_distance(start, end, feature)
            for feature in chromosome_features
        )

        selected = [
            feature for feature in chromosome_features
            if feature_distance(start, end, feature) == minimum
        ]

        relationship = "NEAREST"
        distance = minimum

    return {
        "gene_relationship": relationship,
        "gene_distance_bp": distance,
        "gene_IDs": ",".join(
            sorted({
                feature["gene_ID"]
                for feature in selected
            })
        ),
        "gene_symbols": ",".join(
            sorted({
                feature["gene_symbol"]
                for feature in selected
                if feature["gene_symbol"] not in {"", "."}
            })
        ),
        "feature_types": ",".join(
            sorted({
                feature["feature_type"]
                for feature in selected
            })
        ),
        "locus_types": ",".join(
            sorted({
                feature["locus_type"]
                for feature in selected
            })
        ),
    }


def load_historical_switches(path, label):
    switches = []

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            specific = (
                row["stress_specific_switch"]
                .strip()
                .lower()
            )

            if specific not in {"true", "1", "yes"}:
                continue

            trajectory = row["trajectory"].strip()

            if trajectory == "L>L>H":
                response = "GAIN"
            elif trajectory == "H>H>L":
                response = "LOSS"
            else:
                continue

            switches.append({
                "stress": label,
                "chrom": str(row["chrom"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "response": response,
                "trajectory": trajectory,
            })

    return switches


features = load_features(FEATURES)

heat_switches = load_historical_switches(
    HEAT, "HEAT"
)

cold_switches = load_historical_switches(
    COLD, "COLD"
)

historical_switches = heat_switches + cold_switches


# ------------------------------------------------------------
# Annotate the six recurrent UVC loci
# ------------------------------------------------------------

recurrent_rows = []

with UVC_TABLE.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        overlap_count = int(row["historical_overlap_count"])

        if overlap_count == 0:
            continue

        chrom = str(row["chrom"])
        start = int(row["start"])
        end = int(row["end"])

        annotation = annotate_interval(
            chrom, start, end, features
        )

        output = dict(row)
        output.update(annotation)
        recurrent_rows.append(output)

if len(recurrent_rows) != 6:
    raise RuntimeError(
        f"Expected 6 recurrent UVC loci; "
        f"found {len(recurrent_rows)}"
    )

recurrent_fields = list(recurrent_rows[0])

with RECURRENT_OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=recurrent_fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(recurrent_rows)


# ------------------------------------------------------------
# Audit historical switches inside four 500-bp UVC DMRs
# ------------------------------------------------------------

high_rows = []

with HIGH.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    required = {
        "DMR_chrom", "DMR_start", "DMR_end",
        "DMR_locus", "DMR_direction",
        "treatment_beta", "FDR",
        "strict_locus", "strict_response"
    }

    missing = required - set(reader.fieldnames or [])

    if missing:
        raise RuntimeError(
            f"High-confidence table missing columns: "
            f"{sorted(missing)}"
        )

    for row in reader:
        chrom = str(row["DMR_chrom"])
        start = int(row["DMR_start"])
        end = int(row["DMR_end"])

        hits = [
            switch
            for switch in historical_switches
            if (
                switch["chrom"] == chrom and
                start < switch["end"] and
                end > switch["start"]
            )
        ]

        heat_hits = [
            switch for switch in hits
            if switch["stress"] == "HEAT"
        ]

        cold_hits = [
            switch for switch in hits
            if switch["stress"] == "COLD"
        ]

        annotation = annotate_interval(
            chrom, start, end, features
        )

        output = dict(row)

        output.update({
            "heat_switch_count": len(heat_hits),
            "heat_switch_loci": ",".join(
                f"{hit['chrom']}:{hit['start']}-{hit['end']}"
                for hit in heat_hits
            ),
            "heat_responses": ",".join(
                hit["response"]
                for hit in heat_hits
            ),
            "cold_switch_count": len(cold_hits),
            "cold_switch_loci": ",".join(
                f"{hit['chrom']}:{hit['start']}-{hit['end']}"
                for hit in cold_hits
            ),
            "cold_responses": ",".join(
                hit["response"]
                for hit in cold_hits
            ),
            "historical_switch_count": len(hits),
        })

        output.update(annotation)
        high_rows.append(output)

if len(high_rows) != 4:
    raise RuntimeError(
        f"Expected 4 high-confidence DMRs; "
        f"found {len(high_rows)}"
    )

high_fields = list(high_rows[0])

with HIGH_OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=high_fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(high_rows)


print("===== RECURRENT UVC LOCI =====")
print(f"Records: {len(recurrent_rows)}")

for row in recurrent_rows:
    print(
        row["locus"],
        row["UVC_response"],
        f"heat={row['HEAT_response'] or '-'}",
        f"cold={row['COLD_response'] or '-'}",
        row["gene_relationship"],
        row["gene_IDs"] or "-",
        row["gene_symbols"] or "-",
        sep="\t"
    )

print("\n===== FOUR HIGH-CONFIDENCE UVC DMRs =====")

for row in high_rows:
    print(
        row["DMR_locus"],
        row["DMR_direction"],
        f"heat_switches={row['heat_switch_count']}",
        f"cold_switches={row['cold_switch_count']}",
        row["gene_relationship"],
        row["gene_IDs"] or "-",
        row["gene_symbols"] or "-",
        sep="\t"
    )

print(f"\nWROTE: {RECURRENT_OUT}")
print(f"WROTE: {HIGH_OUT}")
print("RECURRENT AND HIGH-CONFIDENCE LOCUS AUDIT: COMPLETE")
