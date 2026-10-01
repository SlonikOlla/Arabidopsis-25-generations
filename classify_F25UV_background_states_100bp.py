#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
DESIGN = ROOT / "metadata/wgbs_F25C_F25UV_background_adjusted_design.tsv"
WINDOWDIR = ROOT / "results/windows_100bp"
OUTDIR = ROOT / "results/F25UV_background_states_100bp"

CONTEXTS = ["CG", "CHG", "CHH"]

COVERED_CS = {
    "CG": 3,
    "CHG": 3,
    "CHH": 10,
}

STATE_CUTOFFS = {
    "CG": (0.025, 0.10),
    "CHG": (0.025, 0.10),
    "CHH": (0.050, 0.10),
}

design = pd.read_csv(DESIGN, sep="\t", dtype=str)

expected_cells = {
    ("CONTROL", "A"): 3,
    ("UVC", "A"): 2,
    ("CONTROL", "B"): 2,
    ("UVC", "B"): 3,
}

observed = (
    design.groupby(["treatment", "background"])
    .size()
    .to_dict()
)

if observed != expected_cells:
    raise RuntimeError(
        f"Unexpected design cells: {observed}; "
        f"expected {expected_cells}"
    )

OUTDIR.mkdir(parents=True, exist_ok=True)

summary_rows = []
transition_tables = {}

for context in CONTEXTS:

    print(f"===== {context}: READING 10 SAMPLES =====", flush=True)

    series = []

    for row in design.itertuples(index=False):

        path = WINDOWDIR / f"{row.sample}_100bp.tsv"

        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

        data = pd.read_csv(
            path,
            sep="\t",
            usecols=[
                "chrom",
                "start",
                "end",
                "context",
                "covered_Cs",
                "weighted_methylation",
            ],
            dtype={"chrom": str},
        )

        data = data[
            (data["context"] == context) &
            (data["chrom"].isin(["1", "2", "3", "4", "5"])) &
            (data["covered_Cs"] >= COVERED_CS[context])
        ].copy()

        data["locus"] = (
            data["chrom"] + ":" +
            data["start"].astype(str) + "-" +
            data["end"].astype(str)
        )

        values = data.set_index("locus")["weighted_methylation"]
        values.name = row.sample
        series.append(values)

        print(
            row.sample,
            row.treatment,
            row.background,
            len(values),
            flush=True,
        )

    matrix = pd.concat(series, axis=1)

    coordinates = (
        matrix.index.to_series()
        .str.extract(
            r"^(?P<chrom>[^:]+):"
            r"(?P<start>\d+)-(?P<end>\d+)$"
        )
    )

    cell_outputs = {}

    for treatment, background in expected_cells:

        samples = design.loc[
            (design["treatment"] == treatment) &
            (design["background"] == background),
            "sample"
        ].tolist()

        cell = matrix[samples]
        n_expected = len(samples)

        n_covered = cell.notna().sum(axis=1)
        median = cell.median(axis=1, skipna=True)

        low, high = STATE_CUTOFFS[context]

        n_low = (cell < low).sum(axis=1)
        n_high = (cell >= high).sum(axis=1)
        n_buffer = (
            (cell >= low) &
            (cell < high)
        ).sum(axis=1)

        state = pd.Series(
            "AMB",
            index=cell.index,
            dtype=object,
        )

        state[
            (n_covered == n_expected) &
            (n_low == n_expected)
        ] = "L"

        state[
            (n_covered == n_expected) &
            (n_high == n_expected)
        ] = "H"

        cell_name = f"{treatment}_{background}"

        output = pd.DataFrame({
            "chrom": coordinates["chrom"].values,
            "start": coordinates["start"].astype(int).values,
            "end": coordinates["end"].astype(int).values,
            "context": context,
            "cell": cell_name,
            "n_expected": n_expected,
            "n_covered": n_covered.values,
            "n_L": n_low.values,
            "n_H": n_high.values,
            "n_buffer": n_buffer.values,
            "median_methylation": median.values,
            "state": state.values,
        })

        output = output.sort_values(["chrom", "start", "end"])

        outfile = (
            OUTDIR /
            f"{context}_{cell_name}_states_STRICT.tsv.gz"
        )

        output.to_csv(
            outfile,
            sep="\t",
            index=False,
            compression="gzip",
        )

        cell_outputs[cell_name] = output

        counts = output["state"].value_counts()

        for state_name in ["L", "H", "AMB"]:
            summary_rows.append({
                "context": context,
                "cell": cell_name,
                "samples": n_expected,
                "state": state_name,
                "windows": int(counts.get(state_name, 0)),
            })

        print(
            context,
            cell_name,
            dict(counts),
            flush=True,
        )

    for background in ["A", "B"]:

        control = cell_outputs[f"CONTROL_{background}"]
        uvc = cell_outputs[f"UVC_{background}"]

        merged = control.merge(
            uvc,
            on=["chrom", "start", "end", "context"],
            suffixes=("_control", "_UVC"),
            validate="one_to_one",
        )

        confident = (
            merged["state_control"].isin(["L", "H"]) &
            merged["state_UVC"].isin(["L", "H"])
        )

        merged["transition"] = "AMB"
        merged.loc[
            confident,
            "transition"
        ] = (
            merged.loc[confident, "state_control"] +
            merged.loc[confident, "state_UVC"]
        )

        merged["delta_median"] = (
            merged["median_methylation_UVC"] -
            merged["median_methylation_control"]
        )

        outfile = (
            OUTDIR /
            f"{context}_background_{background}_"
            f"CONTROL_to_UVC_transitions.tsv.gz"
        )

        merged.to_csv(
            outfile,
            sep="\t",
            index=False,
            compression="gzip",
        )

        transition_tables[background] = merged

        counts = merged["transition"].value_counts()

        for transition in ["LL", "LH", "HL", "HH", "AMB"]:
            summary_rows.append({
                "context": context,
                "cell": f"TRANSITION_{background}",
                "samples": expected_cells[("CONTROL", background)] +
                           expected_cells[("UVC", background)],
                "state": transition,
                "windows": int(counts.get(transition, 0)),
            })

        print(
            context,
            f"BACKGROUND_{background}",
            dict(counts),
            flush=True,
        )

    a = transition_tables["A"][
        [
            "chrom",
            "start",
            "end",
            "context",
            "transition",
            "delta_median",
        ]
    ].rename(columns={
        "transition": "transition_A",
        "delta_median": "delta_median_A",
    })

    b = transition_tables["B"][
        [
            "chrom",
            "start",
            "end",
            "context",
            "transition",
            "delta_median",
        ]
    ].rename(columns={
        "transition": "transition_B",
        "delta_median": "delta_median_B",
    })

    shared = a.merge(
        b,
        on=["chrom", "start", "end", "context"],
        validate="one_to_one",
    )

    shared["shared_UVC_response"] = "NO"

    shared.loc[
        (shared["transition_A"] == "LH") &
        (shared["transition_B"] == "LH"),
        "shared_UVC_response"
    ] = "GAIN"

    shared.loc[
        (shared["transition_A"] == "HL") &
        (shared["transition_B"] == "HL"),
        "shared_UVC_response"
    ] = "LOSS"

    shared_file = (
        OUTDIR /
        f"{context}_shared_background_concordant_UVC_switches.tsv"
    )

    shared[
        shared["shared_UVC_response"].isin(["GAIN", "LOSS"])
    ].to_csv(
        shared_file,
        sep="\t",
        index=False,
    )

    shared_counts = shared["shared_UVC_response"].value_counts()

    print(
        context,
        "SHARED_RESPONSES",
        dict(shared_counts),
        flush=True,
    )

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUTDIR / "F25UV_background_state_summary.tsv",
    sep="\t",
    index=False,
)

print("F25UV BACKGROUND-STRATIFIED STATES: COMPLETE")
