#!/usr/bin/env bash
#SBATCH --job-name=UVC_WGBS_retry
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --array=1-2%2
#SBATCH --output=logs/UVC/wgbs_production/UVC_WGBS_retry_%A_%a.out
#SBATCH --error=logs/UVC/wgbs_production/UVC_WGBS_retry_%A_%a.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"

samples=(
    F25UV_10_WGBS
    F25UV_11_WGBS
)

sample="${samples[$((SLURM_ARRAY_TASK_ID - 1))]}"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load bismark/0.25.1

echo "========================================"
echo "UVC WGBS RETRY"
echo "SAMPLE: $sample"
echo "TASK:   $SLURM_ARRAY_TASK_ID"
echo "START:  $(date)"
echo "HOST:   $(hostname)"
echo "========================================"

echo "===== SOFTWARE ====="
bismark --version
bowtie2 --version
samtools --version

r1="trimmed/wgbs_all/${sample}_1_val_1.fq.gz"
r2="trimmed/wgbs_all/${sample}_2_val_2.fq.gz"

echo "===== INPUT VALIDATION ====="

for f in "$r1" "$r2"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty input: $f"
        exit 1
    fi
    echo "PASS: $f ($(stat -c '%s' "$f") bytes)"
done

window100="results/windows_100bp/${sample}_100bp.tsv"
window500="results/windows_500bp/${sample}_500bp.tsv"

if [[ -s "$window100" || -s "$window500" ]]; then
    echo "ERROR: final window output already exists for $sample"
    echo "Refusing to overwrite an existing completed analysis."
    exit 1
fi

echo "===== PRODUCTION PIPELINE ====="

bash scripts/process_wgbs_sample.sh "$sample"

if [[ ! -s "$window100" || ! -s "$window500" ]]; then
    echo "ERROR: expected final window files were not produced"
    exit 1
fi

echo "========================================"
echo "UVC WGBS RETRY COMPLETE"
echo "SAMPLE: $sample"
echo "FINISH: $(date)"
echo "========================================"
