#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip

from scipy.stats import fisher_exact

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

UVC_A = (
    ROOT /
    "results/F25UV_background_states_100bp/"
    "CG_background_A_CONTROL_to_UVC_transitions.tsv.gz"
)

UVC_B = (
    ROOT /
    "results/F25UV_background_states_100bp/"
    "CG_background_B_CONTROL_to_UVC_transitions.tsv.gz"
)

OUT = (
    ROOT /
    "results/cross_stress_CG_comparison/"
    "CG_switch_overlap_enrichment.tsv"
)

CONFIDENT = {"LL", "LH", "HL", "HH"}


def locus(row):
    return (
        str(row["chrom"]),
        int(row["start"]),
        int(row["end"]),
    )


def historical_sets(path, label):
    eligible = set()
    switches = {}

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "chrom", "start", "end",
            "trajectory", "stress_specific_switch"
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                f"{label}: missing columns: {sorted(missing)}"
            )

        for row in reader:
            key = locus(row)
            trajectory = row["trajectory"].strip()

            states = trajectory.split(">")

            if (
                len(states) == 3 and
                all(state in {"L", "H"} for state in states)
            ):
                eligible.add(key)

            specific = (
                row["stress_specific_switch"]
                .strip()
                .lower()
            )

            if specific not in {"true", "1", "yes"}:
                continue

            if trajectory == "L>L>H":
                switches[key] = "GAIN"
            elif trajectory == "H>H>L":
                switches[key] = "LOSS"
            else:
                raise RuntimeError(
                    f"{label}: unexpected switch trajectory "
                    f"{trajectory}"
                )

    return eligible, switches


def uvc_transition_sets(path, label):
    eligible = set()
    switches = {}

    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "chrom", "start", "end", "transition"
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                f"{label}: missing columns: {sorted(missing)}"
            )

        for row in reader:
            transition = row["transition"].strip()

            if transition not in CONFIDENT:
                continue

            key = locus(row)
            eligible.add(key)

            if transition == "LH":
                switches[key] = "GAIN"
            elif transition == "HL":
                switches[key] = "LOSS"

    return eligible, switches


heat_eligible, heat_switches = historical_sets(
    HEAT, "HEAT"
)

cold_eligible, cold_switches = historical_sets(
    COLD, "COLD"
)

uvc_a_eligible, uvc_a_switches = uvc_transition_sets(
    UVC_A, "UVC_A"
)

uvc_b_eligible, uvc_b_switches = uvc_transition_sets(
    UVC_B, "UVC_B"
)

uvc_eligible = uvc_a_eligible & uvc_b_eligible

uvc_switches = {
    key: uvc_a_switches[key]
    for key in set(uvc_a_switches) & set(uvc_b_switches)
    if uvc_a_switches[key] == uvc_b_switches[key]
}

if len(heat_switches) != 998:
    raise RuntimeError(
        f"Expected 998 heat switches; found {len(heat_switches)}"
    )

if len(cold_switches) != 386:
    raise RuntimeError(
        f"Expected 386 cold switches; found {len(cold_switches)}"
    )

if len(uvc_switches) != 38:
    raise RuntimeError(
        f"Expected 38 UVC switches; found {len(uvc_switches)}"
    )


def overlap_test(name1, eligible1, switches1,
                 name2, eligible2, switches2):

    universe = eligible1 & eligible2
    set1 = set(switches1) & universe
    set2 = set(switches2) & universe

    both = len(set1 & set2)
    first_only = len(set1 - set2)
    second_only = len(set2 - set1)
    neither = (
        len(universe) -
        both -
        first_only -
        second_only
    )

    odds_ratio, p_upper = fisher_exact(
        [
            [both, first_only],
            [second_only, neither],
        ],
        alternative="greater"
    )

    expected = (
        len(set1) * len(set2) / len(universe)
        if universe else float("nan")
    )

    fold = (
        both / expected
        if expected > 0 else float("nan")
    )

    shared = set1 & set2

    concordant = sum(
        switches1[key] == switches2[key]
        for key in shared
    )

    opposite = both - concordant

    return {
        "comparison": f"{name1}_vs_{name2}",
        "common_eligible_loci": len(universe),
        "first_switches_in_universe": len(set1),
        "second_switches_in_universe": len(set2),
        "observed_shared_switches": both,
        "expected_shared_switches": expected,
        "fold_enrichment": fold,
        "odds_ratio": odds_ratio,
        "Fisher_upper_P": p_upper,
        "direction_concordant": concordant,
        "direction_opposite": opposite,
    }


rows = [
    overlap_test(
        "HEAT", heat_eligible, heat_switches,
        "COLD", cold_eligible, cold_switches
    ),
    overlap_test(
        "HEAT", heat_eligible, heat_switches,
        "UVC", uvc_eligible, uvc_switches
    ),
    overlap_test(
        "COLD", cold_eligible, cold_switches,
        "UVC", uvc_eligible, uvc_switches
    ),
]

fields = [
    "comparison",
    "common_eligible_loci",
    "first_switches_in_universe",
    "second_switches_in_universe",
    "observed_shared_switches",
    "expected_shared_switches",
    "fold_enrichment",
    "odds_ratio",
    "Fisher_upper_P",
    "direction_concordant",
    "direction_opposite",
]

with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()

    for row in rows:
        formatted = row.copy()

        for field in [
            "expected_shared_switches",
            "fold_enrichment",
            "odds_ratio",
            "Fisher_upper_P",
        ]:
            formatted[field] = (
                f"{row[field]:.8g}"
            )

        writer.writerow(formatted)

print("===== ELIGIBLE LOCI =====")
print(f"HEAT: {len(heat_eligible)}")
print(f"COLD: {len(cold_eligible)}")
print(f"UVC A: {len(uvc_a_eligible)}")
print(f"UVC B: {len(uvc_b_eligible)}")
print(f"UVC common A/B: {len(uvc_eligible)}")

print("\n===== OVERLAP ENRICHMENT =====")

for row in rows:
    print(
        f"{row['comparison']}: "
        f"observed={row['observed_shared_switches']} "
        f"expected={row['expected_shared_switches']:.4f} "
        f"fold={row['fold_enrichment']:.3f} "
        f"OR={row['odds_ratio']:.3f} "
        f"P={row['Fisher_upper_P']:.6g} "
        f"concordant={row['direction_concordant']} "
        f"opposite={row['direction_opposite']}"
    )

print(f"\nWROTE: {OUT}")
print("CG SWITCH OVERLAP ENRICHMENT: COMPLETE")
