#!/usr/bin/env bash
#SBATCH --job-name=PlotPCA_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/UVC/primary23/plot_PCA_%j.out
#SBATCH --error=logs/UVC/primary23/plot_PCA_%j.err

set -euo pipefail

cd /scratch/vasilisa/arabidopsis_25gen

module --force purge
module load StdEnv/2023
module load python/3.11
module load scipy-stack

python3 scripts/plot_PRIMARY23_full_PCA.py
