#!/usr/bin/env bash
#SBATCH --job-name=Plot_UVC_WGS
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/UVC/primary23/plot_UVC_WGS_%j.out
#SBATCH --error=logs/UVC/primary23/plot_UVC_WGS_%j.err

set -euo pipefail

cd /scratch/vasilisa/arabidopsis_25gen

module --force purge
module load StdEnv/2023
module load python/3.11
module load scipy-stack

python3 scripts/plot_PRIMARY23_UVC_WGS_summary.py
