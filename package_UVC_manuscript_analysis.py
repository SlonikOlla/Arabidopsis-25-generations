#!/usr/bin/env python3

from pathlib import Path
import fnmatch
import hashlib
import zipfile
import csv

ROOT = Path("/scratch/vasilisa/arabidopsis_25gen")

OUTPUT = (
    ROOT /
    "Arabidopsis25_UVC_analysis_export_20260930.zip"
)

MANIFEST = (
    ROOT /
    "Arabidopsis25_UVC_analysis_export_20260930_manifest.tsv"
)

MAX_FILE_SIZE = 30 * 1024 * 1024

include_patterns = [
    # Metadata and experimental design
    "metadata/wgs_primary23_groups.tsv",
    "metadata/wgbs_all25_qc_status.tsv",
    "metadata/wgbs_primary22_groups.tsv",
    "metadata/wgbs_F25C_F25UV_background_adjusted_design.tsv",

    # WGS PCA results and figures
    "results/wgs_variants/PRIMARY23_PCA/*.variance.tsv",
    "results/wgs_variants/PRIMARY23_PCA/*_annotated.tsv",
    "results/wgs_variants/PRIMARY23_PCA/*_group_summary.tsv",
    "results/wgs_variants/PRIMARY23_PCA/figures/*.png",
    "results/wgs_variants/PRIMARY23_PCA/figures/*.pdf",

    # WGS UVC mutation analysis
    "results/wgs_variants/PRIMARY23_mutation_analysis/*.tsv",
    "results/wgs_variants/PRIMARY23_mutation_analysis/*.csv",
    "results/wgs_variants/PRIMARY23_mutation_analysis/figures/*.png",
    "results/wgs_variants/PRIMARY23_mutation_analysis/figures/*.pdf",

    # Relatedness summaries
    "results/wgs_variants/PRIMARY23_relatedness*/*.tsv",

    # Primary22 WGBS QC and PCA summaries
    "results/wgbs_PRIMARY22_500bp/*summary.tsv",
    "results/wgbs_PRIMARY22_500bp/*PCA_scores.tsv",
    "results/wgbs_PRIMARY22_500bp/*sample_correlations.tsv",
    "results/wgbs_PRIMARY22_500bp/PRIMARY22_sample_global_methylation.tsv",

    # UVC adjusted DMR results, annotations and figures
    "results/wgbs_PRIMARY22_500bp/UVC_background_adjusted/*.tsv",
    "results/wgbs_PRIMARY22_500bp/UVC_background_adjusted/*.bed",
    "results/wgbs_PRIMARY22_500bp/UVC_background_adjusted/figures/*.png",
    "results/wgbs_PRIMARY22_500bp/UVC_background_adjusted/figures/*.pdf",

    # Strict UVC state-switch results
    "results/F25UV_background_states_100bp/*summary.tsv",
    "results/F25UV_background_states_100bp/*shared_background_concordant_UVC_switches.tsv",

    # Cross-stress comparison
    "results/cross_stress_CG_comparison/*.tsv",
    "results/cross_stress_CG_comparison/*.bed",
    "results/cross_stress_CG_comparison/figures/*.png",
    "results/cross_stress_CG_comparison/figures/*.pdf",

    # Relevant scripts
    "scripts/*PRIMARY23*.py",
    "scripts/*PRIMARY23*.sh",
    "scripts/*PRIMARY22*.py",
    "scripts/*PRIMARY22*.sh",
    "scripts/*UVC*.py",
    "scripts/*UVC*.sh",
    "scripts/*F25UV*.py",
    "scripts/*F25UV*.sh",
    "scripts/*heat_cold*.py",
    "scripts/*recurrent*.py",
    "scripts/plot_CG_heat_cold_UVC_comparison.py",
    "scripts/test_CG_switch_overlap_enrichment.py",
    "scripts/annotate_recurrent_and_high_confidence_UVC_loci.py",
    "scripts/describe_recurrent_UVC_stress_genes.py",
]

exclude_patterns = [
    "*.vcf",
    "*.vcf.gz",
    "*.tbi",
    "*.csi",
    "*.bam",
    "*.bai",
    "*.cram",
    "*.pgen",
    "*.pvar",
    "*.psam",
    "*.raw",
    "*.eigenvec.allele",
    "*weighted_methylation.tsv.gz",
    "*complete_genotype_matrix.tsv",
    "*candidate_independent_pileup.tsv",
    "*Araport11_ANN.tsv",
    "*snpEff_summary.csv",
    "Araport11_gene_and_TE_gene_features.bed",
    "*all_tested*.tsv",
    "*background_adjusted_all*.tsv",
]

selected = {}

for path in ROOT.rglob("*"):
    if not path.is_file():
        continue

    relative = path.relative_to(ROOT).as_posix()

    if any(
        fnmatch.fnmatch(relative, pattern)
        for pattern in include_patterns
    ):
        selected[relative] = path

final_files = []

for relative, path in sorted(selected.items()):
    basename = path.name

    if any(
        fnmatch.fnmatch(basename, pattern) or
        fnmatch.fnmatch(relative, pattern)
        for pattern in exclude_patterns
    ):
        continue

    size = path.stat().st_size

    if size > MAX_FILE_SIZE:
        print(
            f"SKIP LARGE: {relative} "
            f"({size / 1024 / 1024:.1f} MB)"
        )
        continue

    final_files.append((relative, path))


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


manifest_rows = []

for relative, path in final_files:
    manifest_rows.append({
        "relative_path": relative,
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    })

with MANIFEST.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "relative_path",
            "size_bytes",
            "sha256",
        ],
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(manifest_rows)

with zipfile.ZipFile(
    OUTPUT,
    mode="w",
    compression=zipfile.ZIP_DEFLATED,
    compresslevel=6
) as archive:

    for relative, path in final_files:
        archive.write(
            path,
            arcname=relative
        )

    archive.write(
        MANIFEST,
        arcname=MANIFEST.name
    )

with zipfile.ZipFile(OUTPUT, "r") as archive:
    bad_file = archive.testzip()

    if bad_file is not None:
        raise RuntimeError(
            f"ZIP integrity failure: {bad_file}"
        )

zip_size = OUTPUT.stat().st_size

print("===== UVC MANUSCRIPT EXPORT =====")
print(f"Files packaged: {len(final_files)}")
print(
    f"ZIP size: {zip_size / 1024 / 1024:.2f} MB"
)
print(f"ZIP integrity: PASS")
print(f"OUTPUT: {OUTPUT}")
print(f"MANIFEST: {MANIFEST}")

print("\n===== FILES BY TOP-LEVEL CATEGORY =====")

categories = {}

for relative, _ in final_files:
    parts = relative.split("/")

    if parts[0] == "results" and len(parts) > 1:
        category = "/".join(parts[:2])
    else:
        category = parts[0]

    categories[category] = categories.get(category, 0) + 1

for category, count in sorted(categories.items()):
    print(f"{category}\t{count}")

print("UVC MANUSCRIPT ANALYSIS EXPORT: COMPLETE")
