#!/usr/bin/env bash
#SBATCH --job-name=WGBS_P22_500
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/UVC/primary22_WGBS_500bp_%j.out
#SBATCH --error=logs/UVC/primary22_WGBS_500bp_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load python/3.11
module load scipy-stack

echo "========================================"
echo "PRIMARY22 500-BP WGBS MATRIX BUILD"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

python3 scripts/build_PRIMARY22_WGBS_500bp_matrices.py

echo "FINISH: $(date)"
