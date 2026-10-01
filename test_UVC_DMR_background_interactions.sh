#!/usr/bin/env bash
#SBATCH --job-name=UVC_DMR_INT
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/wgbs_background/UVC_DMR_interaction_%j.out
#SBATCH --error=logs/UVC/wgbs_background/UVC_DMR_interaction_%j.err

set -euo pipefail

cd /scratch/vasilisa/arabidopsis_25gen

module --force purge
module load StdEnv/2023
module load python/3.11
module load scipy-stack

python3 scripts/test_UVC_DMR_background_interactions.py
