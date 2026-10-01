#!/usr/bin/env python3

from pathlib import Path
import csv
import math
import sys

import numpy as np
import pandas as pd

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
META = ROOT / "metadata/wgbs_primary22_groups.tsv"
WINDOWDIR = ROOT / "results/windows_500bp"
OUTDIR = ROOT / "results/wgbs_PRIMARY22_500bp"

CONTEXTS = ("CG", "CHG", "CHH")
THRESHOLDS = {"CG": 3, "CHG": 3, "CHH": 10}
NUCLEAR = {"1", "2", "3", "4", "5"}

OUTDIR.mkdir(parents=True, exist_ok=True)

metadata = pd.read_csv(META, sep="\t", dtype=str)

if list(metadata.columns) != ["sample", "group"]:
    raise RuntimeError(
        f"Unexpected metadata columns: {list(metadata.columns)}"
    )

if len(metadata) != 22:
    raise RuntimeError(
        f"Expected 22 metadata records but found {len(metadata)}"
    )

if metadata["sample"].duplicated().any():
    raise RuntimeError("Duplicate samples found in metadata")

samples = metadata["sample"].tolist()
retain_n = math.ceil(0.80 * len(samples))

print(f"SAMPLES: {len(samples)}", flush=True)
print(f"80% RETENTION COUNT: {retain_n}", flush=True)

global_rows = []
context_series = {context: [] for context in CONTEXTS}
input_rows = {context: {} for context in CONTEXTS}
passing_rows = {context: {} for context in CONTEXTS}

for sample in samples:
    path = WINDOWDIR / f"{sample}_500bp.tsv"

    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)

    print(f"READING: {sample}", flush=True)

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=[
            "chrom", "start", "end", "context",
            "meth_reads", "total_reads",
            "covered_Cs", "weighted_methylation"
        ],
        dtype={
            "chrom": str,
            "start": np.int32,
            "end": np.int32,
            "context": str,
            "meth_reads": np.int64,
            "total_reads": np.int64,
            "covered_Cs": np.int32,
            "weighted_methylation": np.float32
        }
    )

    df = df[df["chrom"].isin(NUCLEAR)].copy()

    for context in CONTEXTS:
        sub = df[df["context"] == context].copy()
        input_rows[context][sample] = len(sub)

        meth_sum = int(sub["meth_reads"].sum())
        total_sum = int(sub["total_reads"].sum())
        covered_sum = int(sub["covered_Cs"].sum())

        weighted = (
            meth_sum / total_sum
            if total_sum > 0
            else np.nan
        )

        global_rows.append({
            "sample": sample,
            "group": metadata.loc[
                metadata["sample"] == sample, "group"
            ].iloc[0],
            "context": context,
            "meth_reads": meth_sum,
            "total_reads": total_sum,
            "covered_Cs_sum": covered_sum,
            "global_weighted_methylation": weighted
        })

        sub = sub[
            sub["covered_Cs"] >= THRESHOLDS[context]
        ].copy()

        passing_rows[context][sample] = len(sub)

        sub["locus"] = (
            sub["chrom"] + ":" +
            sub["start"].astype(str) + "-" +
            sub["end"].astype(str)
        )

        if sub["locus"].duplicated().any():
            raise RuntimeError(
                f"Duplicate {context} loci found in {sample}"
            )

        series = sub.set_index("locus")["weighted_methylation"]
        series.name = sample
        context_series[context].append(series)

    del df

global_df = pd.DataFrame(global_rows)
global_df.to_csv(
    OUTDIR / "PRIMARY22_sample_global_methylation.tsv",
    sep="\t",
    index=False,
    float_format="%.8f"
)

summary_rows = []

for context in CONTEXTS:
    print(f"BUILDING MATRIX: {context}", flush=True)

    matrix = pd.concat(
        context_series[context],
        axis=1,
        join="outer"
    )

    observed = matrix.notna().sum(axis=1)
    matrix = matrix.loc[observed >= retain_n].copy()
    observed = observed.loc[matrix.index]

    matrix.insert(0, "samples_observed", observed.astype(int))

    matrix.to_csv(
        OUTDIR /
        f"PRIMARY22_{context}_500bp_weighted_methylation.tsv.gz",
        sep="\t",
        compression="gzip",
        float_format="%.6f"
    )

    scientific_matrix = matrix.drop(
        columns=["samples_observed"]
    )

    correlations = scientific_matrix.corr(
        method="pearson",
        min_periods=100
    )

    correlations.to_csv(
        OUTDIR / f"PRIMARY22_{context}_sample_correlations.tsv",
        sep="\t",
        float_format="%.8f"
    )

    pca_matrix = scientific_matrix.apply(
        lambda row: row.fillna(row.median()),
        axis=1
    )

    variable = pca_matrix.var(axis=1) > 0
    pca_matrix = pca_matrix.loc[variable]

    X = pca_matrix.T.to_numpy(dtype=np.float64, copy=True)
    X -= X.mean(axis=0)

    U, singular_values, _ = np.linalg.svd(
        X,
        full_matrices=False
    )

    scores = U * singular_values
    variance = singular_values ** 2
    explained = variance / variance.sum()

    pc_count = min(10, scores.shape[1])

    pca = pd.DataFrame(
        scores[:, :pc_count],
        index=pca_matrix.columns,
        columns=[
            f"PC{i}" for i in range(1, pc_count + 1)
        ]
    )

    pca.insert(
        0,
        "group",
        metadata.set_index("sample")
        .loc[pca.index, "group"]
    )

    for i in range(pc_count):
        pca[f"PC{i+1}_variance_pct"] = explained[i] * 100

    pca.to_csv(
        OUTDIR / f"PRIMARY22_{context}_PCA_scores.tsv",
        sep="\t",
        index_label="sample",
        float_format="%.8f"
    )

    for sample in samples:
        summary_rows.append({
            "context": context,
            "sample": sample,
            "input_nuclear_windows":
                input_rows[context][sample],
            "windows_passing_covered_C_threshold":
                passing_rows[context][sample],
            "covered_C_threshold":
                THRESHOLDS[context],
            "retained_matrix_loci":
                len(scientific_matrix),
            "variable_PCA_loci":
                len(pca_matrix),
            "minimum_samples_required":
                retain_n
        })

    print(
        f"{context}: retained={len(scientific_matrix)} "
        f"variable={len(pca_matrix)} "
        f"PC1={explained[0]*100:.2f}% "
        f"PC2={explained[1]*100:.2f}%",
        flush=True
    )

summary = pd.DataFrame(summary_rows)
summary.to_csv(
    OUTDIR / "PRIMARY22_500bp_matrix_summary.tsv",
    sep="\t",
    index=False
)

print("PRIMARY22 500-BP WGBS MATRICES: COMPLETE", flush=True)
