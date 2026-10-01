#!/usr/bin/env python3

from pathlib import Path
import itertools

import numpy as np
import pandas as pd
from scipy.stats import (
    kruskal,
    mannwhitneyu,
    ttest_ind
)

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")
INDIR = ROOT / "results/wgbs_PRIMARY22_500bp"
META = ROOT / "metadata/wgbs_primary22_groups.tsv"

GROUP_ORDER = ["F2C", "F25C", "F25H", "F25Cd", "F25UV"]
CONTEXTS = ["CG", "CHG", "CHH"]

metadata = pd.read_csv(META, sep="\t")
metadata["group"] = pd.Categorical(
    metadata["group"],
    categories=GROUP_ORDER,
    ordered=True
)

global_df = pd.read_csv(
    INDIR / "PRIMARY22_sample_global_methylation.tsv",
    sep="\t"
)

if len(global_df) != 66:
    raise RuntimeError(
        f"Expected 66 sample-context records, found {len(global_df)}"
    )

if global_df["sample"].nunique() != 22:
    raise RuntimeError("Expected 22 unique samples")

group_rows = []
test_rows = []

for context in CONTEXTS:
    context_df = global_df[
        global_df["context"] == context
    ].copy()

    for group in GROUP_ORDER:
        values = context_df.loc[
            context_df["group"] == group,
            "global_weighted_methylation"
        ].dropna().to_numpy()

        n = len(values)
        sd = np.std(values, ddof=1) if n > 1 else np.nan
        se = sd / np.sqrt(n) if n > 1 else np.nan

        group_rows.append({
            "context": context,
            "group": group,
            "n": n,
            "mean": np.mean(values),
            "SD": sd,
            "SE": se,
            "minimum": np.min(values),
            "maximum": np.max(values)
        })

    arrays = [
        context_df.loc[
            context_df["group"] == group,
            "global_weighted_methylation"
        ].dropna().to_numpy()
        for group in GROUP_ORDER
    ]

    kw = kruskal(*arrays)

    uvc = context_df.loc[
        context_df["group"] == "F25UV",
        "global_weighted_methylation"
    ].dropna().to_numpy()

    control = context_df.loc[
        context_df["group"] == "F25C",
        "global_weighted_methylation"
    ].dropna().to_numpy()

    welch = ttest_ind(
        uvc,
        control,
        equal_var=False
    )

    mann_whitney = mannwhitneyu(
        uvc,
        control,
        alternative="two-sided"
    )

    test_rows.append({
        "context": context,
        "comparison": "all_five_groups",
        "test": "Kruskal_Wallis",
        "statistic": kw.statistic,
        "P_value": kw.pvalue,
        "mean_difference": np.nan,
        "details": "F2C,F25C,F25H,F25Cd,F25UV"
    })

    test_rows.append({
        "context": context,
        "comparison": "F25UV_vs_F25C",
        "test": "Welch_t_test",
        "statistic": welch.statistic,
        "P_value": welch.pvalue,
        "mean_difference": np.mean(uvc) - np.mean(control),
        "details": "positive_difference=F25UV_higher"
    })

    test_rows.append({
        "context": context,
        "comparison": "F25UV_vs_F25C",
        "test": "Mann_Whitney_U",
        "statistic": mann_whitney.statistic,
        "P_value": mann_whitney.pvalue,
        "mean_difference": np.mean(uvc) - np.mean(control),
        "details": "two_sided"
    })

group_summary = pd.DataFrame(group_rows)
group_summary.to_csv(
    INDIR / "PRIMARY22_global_methylation_group_summary.tsv",
    sep="\t",
    index=False,
    float_format="%.8f"
)

tests = pd.DataFrame(test_rows)
tests.to_csv(
    INDIR / "PRIMARY22_global_methylation_tests.tsv",
    sep="\t",
    index=False,
    float_format="%.8g"
)

pca_tables = []
pca_group_rows = []

for context in CONTEXTS:
    pca = pd.read_csv(
        INDIR / f"PRIMARY22_{context}_PCA_scores.tsv",
        sep="\t"
    )

    if len(pca) != 22:
        raise RuntimeError(
            f"{context}: expected 22 PCA samples, found {len(pca)}"
        )

    pca.insert(0, "context", context)
    pca_tables.append(pca)

    for group in GROUP_ORDER:
        sub = pca[pca["group"] == group]

        row = {
            "context": context,
            "group": group,
            "n": len(sub)
        }

        for pc in ["PC1", "PC2", "PC3"]:
            row[f"{pc}_mean"] = sub[pc].mean()
            row[f"{pc}_SD"] = sub[pc].std(ddof=1)

        pca_group_rows.append(row)

pca_all = pd.concat(pca_tables, ignore_index=True)
pca_all.to_csv(
    INDIR / "PRIMARY22_PCA_annotated_coordinates.tsv",
    sep="\t",
    index=False,
    float_format="%.8f"
)

pca_group = pd.DataFrame(pca_group_rows)
pca_group.to_csv(
    INDIR / "PRIMARY22_PCA_group_summary.tsv",
    sep="\t",
    index=False,
    float_format="%.8f"
)

pair_rows = []

sample_to_group = dict(
    zip(metadata["sample"], metadata["group"].astype(str))
)

for context in CONTEXTS:
    corr = pd.read_csv(
        INDIR / f"PRIMARY22_{context}_sample_correlations.tsv",
        sep="\t",
        index_col=0
    )

    if corr.shape != (22, 22):
        raise RuntimeError(
            f"{context}: correlation matrix shape is {corr.shape}"
        )

    for sample1, sample2 in itertools.combinations(corr.columns, 2):
        group1 = sample_to_group[sample1]
        group2 = sample_to_group[sample2]

        if group1 == group2:
            relationship = "within_group"
            group_pair = group1
        else:
            relationship = "between_group"
            group_pair = "|".join(sorted([group1, group2]))

        pair_rows.append({
            "context": context,
            "sample1": sample1,
            "group1": group1,
            "sample2": sample2,
            "group2": group2,
            "relationship": relationship,
            "group_pair": group_pair,
            "correlation": corr.loc[sample1, sample2]
        })

pairs = pd.DataFrame(pair_rows)
pairs.to_csv(
    INDIR / "PRIMARY22_pairwise_sample_correlations.tsv",
    sep="\t",
    index=False,
    float_format="%.8f"
)

corr_summary = (
    pairs.groupby(
        ["context", "relationship", "group_pair"],
        as_index=False
    )["correlation"]
    .agg(
        n="count",
        mean="mean",
        SD="std",
        minimum="min",
        maximum="max"
    )
)

corr_summary.to_csv(
    INDIR / "PRIMARY22_correlation_group_summary.tsv",
    sep="\t",
    index=False,
    float_format="%.8f"
)

print("PRIMARY22 WGBS QC SUMMARY: COMPLETE")
print(f"Global records: {len(global_df)}")
print(f"PCA records: {len(pca_all)}")
print(f"Correlation pairs: {len(pairs)}")
