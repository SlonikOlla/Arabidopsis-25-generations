#!/usr/bin/env bash
#SBATCH --job-name=UVC_BG_500
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=logs/UVC/wgbs_background/UVC_BG_500_%j.out
#SBATCH --error=logs/UVC/wgbs_background/UVC_BG_500_%j.err

set -euo pipefail

cd /scratch/vasilisa/arabidopsis_25gen

module --force purge
module load StdEnv/2023
module load python/3.11
module load scipy-stack

python3 scripts/test_UVC_background_adjusted_500bp.py
