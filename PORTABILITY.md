
# Portability and reproducibility

## Python environment

Use `environment.yml` or `requirements.txt` for Python dependencies.

## External software used by HPC workflows

- BCFtools 1.22
- BEDTools 2.31.0
- Bismark 0.25.1
- Bowtie2 2.5.x
- GATK
- PLINK 2.0.0-a.6.32
- SnpEff 5.2 with a custom Araport11/TAIR10 database
- Slurm

## Before running

Search scripts for `/scratch/vasilisa/arabidopsis_25gen`, `#SBATCH`, and `module load`. Replace cluster paths, account names, reference paths, software modules, and requested resources as appropriate. The repository does not contain input data, manuscripts, figures, or results.

## Suggested hardening before public release

1. Replace absolute paths with command-line arguments or a YAML configuration.
2. Add a small synthetic test dataset.
3. Add continuous integration for Python syntax and shell linting.
4. Record exact GATK and reference-assembly versions.
5. Link the immutable data archive and article DOI.
