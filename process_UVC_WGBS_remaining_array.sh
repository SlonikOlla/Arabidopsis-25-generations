#!/usr/bin/env bash
#SBATCH --job-name=UVC_WGBS_prod
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --array=1-4%2
#SBATCH --output=logs/UVC/wgbs_production/UVC_WGBS_prod_%A_%a.out
#SBATCH --error=logs/UVC/wgbs_production/UVC_WGBS_prod_%A_%a.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"

samples=(
    F25UV_8_WGBS
    F25UV_10_WGBS
    F25UV_11_WGBS
    F25UV_13_WGBS
)

sample="${samples[$((SLURM_ARRAY_TASK_ID - 1))]}"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load bismark/0.25.1

echo "========================================"
echo "UVC WGBS PRODUCTION"
echo "SAMPLE: $sample"
echo "TASK:   ${SLURM_ARRAY_TASK_ID}"
echo "START:  $(date)"
echo "HOST:   $(hostname)"
echo "========================================"

echo "===== SOFTWARE ====="
bismark --version 2>&1 | grep -m1 -E 'Bismark Version|version'
bowtie2 --version 2>&1 | head -1
samtools --version | head -1

R1="trimmed/wgbs_all/${sample}_1_val_1.fq.gz"
R2="trimmed/wgbs_all/${sample}_2_val_2.fq.gz"

echo "===== INPUT VALIDATION ====="

for f in "$R1" "$R2"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty input: $f"
        exit 1
    fi
done

echo "R1: $(stat -c '%s bytes' "$R1")"
echo "R2: $(stat -c '%s bytes' "$R2")"
echo "INPUTS: PASS"

bash scripts/process_wgbs_sample.sh "$sample"

echo "========================================"
echo "UVC WGBS ARRAY TASK COMPLETE"
echo "SAMPLE: $sample"
echo "FINISH: $(date)"
echo "========================================"
