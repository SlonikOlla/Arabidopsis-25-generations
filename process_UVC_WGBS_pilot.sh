#!/usr/bin/env bash
#SBATCH --job-name=UVC_WGBS_pilot
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=logs/UVC/wgbs_production/UVC_WGBS_pilot_%j.out
#SBATCH --error=logs/UVC/wgbs_production/UVC_WGBS_pilot_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
SAMPLE="F25UV_7_WGBS"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load bismark/0.25.1

echo "========================================"
echo "UVC WGBS PRODUCTION PILOT"
echo "SAMPLE: $SAMPLE"
echo "START:  $(date)"
echo "HOST:   $(hostname)"
echo "========================================"

echo "===== SOFTWARE ====="
bismark --version 2>&1 | grep -m1 -E 'Bismark Version|version'
bowtie2 --version 2>&1 | head -1
samtools --version | head -1

R1="trimmed/wgbs_all/${SAMPLE}_1_val_1.fq.gz"
R2="trimmed/wgbs_all/${SAMPLE}_2_val_2.fq.gz"

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

echo "===== PRODUCTION PIPELINE ====="

bash scripts/process_wgbs_sample.sh "$SAMPLE"

echo "========================================"
echo "UVC WGBS PILOT COMPLETE"
echo "SAMPLE: $SAMPLE"
echo "FINISH: $(date)"
echo "========================================"
