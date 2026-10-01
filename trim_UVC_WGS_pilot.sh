#!/usr/bin/env bash
#SBATCH --job-name=UVC_WGS_trim7
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/wgs_trimmed/UVC_WGS_trim7_%j.out
#SBATCH --error=logs/UVC/wgs_trimmed/UVC_WGS_trim7_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
SAMPLE="F25UV_7_WGS"
TRIM="${ROOT}/tools/trim_galore-2.3.0/trim_galore"
R1="${ROOT}/raw/UVC_normalized/${SAMPLE}_1.fastq.gz"
R2="${ROOT}/raw/UVC_normalized/${SAMPLE}_2.fastq.gz"
OUTDIR="${ROOT}/trimmed/UVC_WGS"

cd "$ROOT"

echo "========================================"
echo "UVC WGS TRIMMING PILOT"
echo "SAMPLE: $SAMPLE"
echo "START:  $(date)"
echo "HOST:   $(hostname)"
echo "========================================"

test -x "$TRIM"
test -s "$R1"
test -s "$R2"

echo "===== SOFTWARE ====="
"$TRIM" --version

echo "===== INPUTS ====="
printf "R1: %s bytes\n" "$(stat -Lc '%s' "$R1")"
printf "R2: %s bytes\n" "$(stat -Lc '%s' "$R2")"

"$TRIM" \
  --paired \
  --cores 4 \
  --quality 20 \
  --length 30 \
  --gzip \
  --output_dir "$OUTDIR" \
  "$R1" "$R2"

TRIM_R1="${OUTDIR}/${SAMPLE}_1_val_1.fq.gz"
TRIM_R2="${OUTDIR}/${SAMPLE}_2_val_2.fq.gz"

test -s "$TRIM_R1"
test -s "$TRIM_R2"

gzip -t "$TRIM_R1"
gzip -t "$TRIM_R2"

echo "===== OUTPUTS VALIDATED ====="
printf "R1: %s bytes\n" "$(stat -c '%s' "$TRIM_R1")"
printf "R2: %s bytes\n" "$(stat -c '%s' "$TRIM_R2")"

echo "RAW INPUTS RETAINED"
echo "UVC WGS TRIMMING PILOT: COMPLETE"
echo "FINISH: $(date)"
