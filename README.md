
# Arabidopsis 25-generation stress inheritance: UVC analysis code

[![Release](https://img.shields.io/badge/release-v1.0.0-blue)](RELEASE_NOTES.md)
[![Status](https://img.shields.io/badge/status-release%20candidate-orange)](PRE_PUBLICATION_CHECKLIST.md)

Analysis-code package for a multigenerational *Arabidopsis thaliana* study comparing control, heat, cold, and UVC lineages with whole-genome sequencing (WGS) and whole-genome bisulfite sequencing (WGBS).

**This package intentionally excludes manuscripts, figures, and result tables.** It contains analysis scripts, minimal sample-design metadata, and release/reproducibility documentation.

## Layout

| Path | Contents |
| --- | --- |
| `scripts/` | WGS, WGBS, annotation, statistics, and plotting source code |
| `metadata/` | Sample groups, WGBS QC status, and UVC background design |
| `docs/` | Workflow and portability notes |
| `checksums/` | SHA-256 checksums |

## Quick validation

```bash
conda env create -f environment.yml
conda activate arabidopsis25-uvc
python -m py_compile scripts/*.py
sha256sum -c checksums/SHA256SUMS
```

## Running analyses

The scripts were developed for a Slurm HPC environment and some contain cluster-specific absolute paths, accounts, modules, and resource directives. Review and parameterize those values before running. Input data and result directories must be supplied separately from the associated data archive.

## Citation and license

Update [`CITATION.cff`](CITATION.cff) with the final author list, ORCIDs, repository URL, and DOI. No reuse license has yet been granted; see [`LICENSE_PENDING.md`](LICENSE_PENDING.md).
