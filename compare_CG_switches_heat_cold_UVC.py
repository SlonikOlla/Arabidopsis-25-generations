#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")

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

UVC = (
    ROOT /
    "results/F25UV_background_states_100bp/"
    "CG_shared_background_concordant_UVC_switches.tsv"
)

OUTDIR = ROOT / "results/cross_stress_CG_comparison"
OUTDIR.mkdir(parents=True, exist_ok=True)

MASTER = OUTDIR / "CG_heat_cold_UVC_switch_master.tsv"
UVC_COMPARISON = OUTDIR / "CG_UVC_switches_heat_cold_comparison.tsv"
PAIRWISE = OUTDIR / "CG_heat_cold_UVC_pairwise_overlap_summary.tsv"
PATTERNS = OUTDIR / "CG_heat_cold_UVC_overlap_patterns.tsv"


def key_from_row(row):
    return (
        str(row["chrom"]),
        int(row["start"]),
        int(row["end"]),
    )


def read_historical(path, stress):
    records = {}

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "chrom", "start", "end",
            "trajectory", "stress_specific_switch"
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                f"{stress}: missing columns: {sorted(missing)}"
            )

        for row in reader:
            specific = row["stress_specific_switch"].strip().lower()

            if specific not in {"true", "1", "yes"}:
                continue

            trajectory = row["trajectory"].strip()

            if trajectory == "L>L>H":
                response = "GAIN"
            elif trajectory == "H>H>L":
                response = "LOSS"
            else:
                raise RuntimeError(
                    f"{stress}: unexpected stress-specific "
                    f"trajectory: {trajectory}"
                )

            key = key_from_row(row)

            if key in records:
                raise RuntimeError(
                    f"{stress}: duplicate locus: {key}"
                )

            records[key] = response

    return records


def read_uvc(path):
    records = {}
    details = {}

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "chrom", "start", "end",
            "shared_UVC_response",
            "delta_median_A", "delta_median_B"
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                f"UVC: missing columns: {sorted(missing)}"
            )

        for row in reader:
            key = key_from_row(row)
            response = row["shared_UVC_response"].strip()

            if response not in {"GAIN", "LOSS"}:
                raise RuntimeError(
                    f"UVC: invalid response at {key}: {response}"
                )

            if key in records:
                raise RuntimeError(
                    f"UVC: duplicate locus: {key}"
                )

            records[key] = response
            details[key] = {
                "delta_A": float(row["delta_median_A"]),
                "delta_B": float(row["delta_median_B"]),
            }

    return records, details


def relationship(first, second):
    if first == "" or second == "":
        return ""

    if first == second:
        return "CONCORDANT"

    return "OPPOSITE"


heat = read_historical(HEAT, "HEAT")
cold = read_historical(COLD, "COLD")
uvc, uvc_details = read_uvc(UVC)

expected = {
    "HEAT": 998,
    "COLD": 386,
    "UVC": 38,
}

observed = {
    "HEAT": len(heat),
    "COLD": len(cold),
    "UVC": len(uvc),
}

for stress, expected_count in expected.items():
    if observed[stress] != expected_count:
        raise RuntimeError(
            f"{stress}: expected {expected_count} switches, "
            f"found {observed[stress]}"
        )

all_loci = sorted(
    set(heat) | set(cold) | set(uvc),
    key=lambda x: (int(x[0]), x[1], x[2])
)

master_fields = [
    "chrom", "start", "end", "locus",
    "HEAT_response", "COLD_response", "UVC_response",
    "HEAT_COLD_relationship",
    "HEAT_UVC_relationship",
    "COLD_UVC_relationship",
    "stress_count",
    "response_pattern",
]

with MASTER.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=master_fields,
        delimiter="\t"
    )
    writer.writeheader()

    for key in all_loci:
        chrom, start, end = key

        h = heat.get(key, "")
        c = cold.get(key, "")
        u = uvc.get(key, "")

        present = [
            stress
            for stress, response in [
                ("HEAT", h),
                ("COLD", c),
                ("UVC", u),
            ]
            if response
        ]

        pattern = ";".join(
            f"{stress}:{response}"
            for stress, response in [
                ("HEAT", h),
                ("COLD", c),
                ("UVC", u),
            ]
            if response
        )

        writer.writerow({
            "chrom": chrom,
            "start": start,
            "end": end,
            "locus": f"{chrom}:{start}-{end}",
            "HEAT_response": h,
            "COLD_response": c,
            "UVC_response": u,
            "HEAT_COLD_relationship": relationship(h, c),
            "HEAT_UVC_relationship": relationship(h, u),
            "COLD_UVC_relationship": relationship(c, u),
            "stress_count": len(present),
            "response_pattern": pattern,
        })


uvc_fields = [
    "chrom", "start", "end", "locus",
    "UVC_response", "UVC_delta_A", "UVC_delta_B",
    "HEAT_response", "HEAT_UVC_relationship",
    "COLD_response", "COLD_UVC_relationship",
    "historical_overlap_count",
]

with UVC_COMPARISON.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=uvc_fields,
        delimiter="\t"
    )
    writer.writeheader()

    for key in sorted(
        uvc,
        key=lambda x: (int(x[0]), x[1], x[2])
    ):
        chrom, start, end = key
        h = heat.get(key, "")
        c = cold.get(key, "")
        u = uvc[key]

        writer.writerow({
            "chrom": chrom,
            "start": start,
            "end": end,
            "locus": f"{chrom}:{start}-{end}",
            "UVC_response": u,
            "UVC_delta_A": uvc_details[key]["delta_A"],
            "UVC_delta_B": uvc_details[key]["delta_B"],
            "HEAT_response": h,
            "HEAT_UVC_relationship": relationship(h, u),
            "COLD_response": c,
            "COLD_UVC_relationship": relationship(c, u),
            "historical_overlap_count":
                int(bool(h)) + int(bool(c)),
        })


def pair_stats(first_name, first, second_name, second):
    shared = set(first) & set(second)
    concordant = sum(
        first[key] == second[key]
        for key in shared
    )
    opposite = len(shared) - concordant

    return {
        "comparison": f"{first_name}_vs_{second_name}",
        "first_total": len(first),
        "second_total": len(second),
        "shared_loci": len(shared),
        "concordant": concordant,
        "opposite": opposite,
        "first_overlap_percent":
            100 * len(shared) / len(first),
        "second_overlap_percent":
            100 * len(shared) / len(second),
    }


pair_rows = [
    pair_stats("HEAT", heat, "COLD", cold),
    pair_stats("HEAT", heat, "UVC", uvc),
    pair_stats("COLD", cold, "UVC", uvc),
]

with PAIRWISE.open("w", newline="") as handle:
    fields = list(pair_rows[0])
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()

    for row in pair_rows:
        row = row.copy()
        row["first_overlap_percent"] = (
            f"{row['first_overlap_percent']:.4f}"
        )
        row["second_overlap_percent"] = (
            f"{row['second_overlap_percent']:.4f}"
        )
        writer.writerow(row)


pattern_counts = Counter()

for key in all_loci:
    h = heat.get(key, "")
    c = cold.get(key, "")
    u = uvc.get(key, "")

    pattern = ";".join(
        f"{stress}:{response}"
        for stress, response in [
            ("HEAT", h),
            ("COLD", c),
            ("UVC", u),
        ]
        if response
    )

    pattern_counts[pattern] += 1

with PATTERNS.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["response_pattern", "loci"])

    for pattern, count in sorted(
        pattern_counts.items(),
        key=lambda item: (-item[1], item[0])
    ):
        writer.writerow([pattern, count])


three_way = set(heat) & set(cold) & set(uvc)
uvc_heat = set(uvc) & set(heat)
uvc_cold = set(uvc) & set(cold)
uvc_any = set(uvc) & (set(heat) | set(cold))

print("===== INPUT SWITCH SETS =====")
print(f"HEAT: {len(heat)}")
print(f"  GAIN: {sum(x == 'GAIN' for x in heat.values())}")
print(f"  LOSS: {sum(x == 'LOSS' for x in heat.values())}")
print(f"COLD: {len(cold)}")
print(f"  GAIN: {sum(x == 'GAIN' for x in cold.values())}")
print(f"  LOSS: {sum(x == 'LOSS' for x in cold.values())}")
print(f"UVC: {len(uvc)}")
print(f"  GAIN: {sum(x == 'GAIN' for x in uvc.values())}")
print(f"  LOSS: {sum(x == 'LOSS' for x in uvc.values())}")

print("\n===== PAIRWISE OVERLAPS =====")
for row in pair_rows:
    print(
        f"{row['comparison']}: "
        f"shared={row['shared_loci']} "
        f"concordant={row['concordant']} "
        f"opposite={row['opposite']}"
    )

print("\n===== UVC HISTORICAL OVERLAP =====")
print(f"UVC overlapping HEAT: {len(uvc_heat)}")
print(f"UVC overlapping COLD: {len(uvc_cold)}")
print(f"UVC overlapping either: {len(uvc_any)}")
print(f"UVC unique: {len(uvc) - len(uvc_any)}")
print(f"Three-way overlap: {len(three_way)}")

print("\n===== OUTPUTS =====")
for path in [MASTER, UVC_COMPARISON, PAIRWISE, PATTERNS]:
    print(path)

print("CROSS-STRESS CG SWITCH COMPARISON: COMPLETE")
