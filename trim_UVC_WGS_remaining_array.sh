#!/usr/bin/env bash
#SBATCH --job-name=UVC_WGS_trim
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-4%2
#SBATCH --output=logs/UVC/wgs_trimmed/UVC_WGS_trim_%A_%a.out
#SBATCH --error=logs/UVC/wgs_trimmed/UVC_WGS_trim_%A_%a.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
TRIM="${ROOT}/tools/trim_galore-2.3.0/trim_galore"
OUTDIR="${ROOT}/trimmed/UVC_WGS"

samples=(
    F25UV_8_WGS
    F25UV_10_WGS
    F25UV_11_WGS
    F25UV_13_WGS
)

SAMPLE="${samples[$((SLURM_ARRAY_TASK_ID - 1))]}"
R1="${ROOT}/raw/UVC_normalized/${SAMPLE}_1.fastq.gz"
R2="${ROOT}/raw/UVC_normalized/${SAMPLE}_2.fastq.gz"

cd "$ROOT"

echo "========================================"
echo "UVC WGS TRIMMING"
echo "SAMPLE: $SAMPLE"
echo "TASK:   $SLURM_ARRAY_TASK_ID"
echo "START:  $(date)"
echo "HOST:   $(hostname)"
echo "========================================"

test -x "$TRIM"
test -s "$R1"
test -s "$R2"

printf "R1 INPUT: %s bytes\n" "$(stat -Lc '%s' "$R1")"
printf "R2 INPUT: %s bytes\n" "$(stat -Lc '%s' "$R2")"

"$TRIM" \
  --paired \
  --cores 4 \
  --quality 20 \
  --length 30 \
  --output_dir "$OUTDIR" \
  "$R1" "$R2"

TRIM_R1="${OUTDIR}/${SAMPLE}_1_val_1.fq.gz"
TRIM_R2="${OUTDIR}/${SAMPLE}_2_val_2.fq.gz"

test -s "$TRIM_R1"
test -s "$TRIM_R2"

gzip -t "$TRIM_R1"
gzip -t "$TRIM_R2"

echo "===== OUTPUTS VALIDATED ====="
printf "R1 OUTPUT: %s bytes\n" "$(stat -c '%s' "$TRIM_R1")"
printf "R2 OUTPUT: %s bytes\n" "$(stat -c '%s' "$TRIM_R2")"

echo "RAW INPUTS RETAINED"
echo "COMPLETE: $SAMPLE"
echo "FINISH: $(date)"
