#!/usr/bin/env bash
#SBATCH --job-name=UVC_WGS_FastQC
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/UVC/wgs_trimmed_fastqc/UVC_WGS_fastqc_%j.out
#SBATCH --error=logs/UVC/wgs_trimmed_fastqc/UVC_WGS_fastqc_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
INPUT="${ROOT}/trimmed/UVC_WGS"
OUTPUT="${ROOT}/qc/UVC_WGS_trimmed_fastqc"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load fastqc/0.12.1

export JAVA_TOOL_OPTIONS="-Xmx2g"

echo "========================================"
echo "UVC WGS POST-TRIM FASTQC"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

mapfile -t files < <(
    find "$INPUT" -maxdepth 1 -type f \
      -name 'F25UV_*_WGS_*_val_*.fq.gz' |
    sort
)

echo "INPUT FILES: ${#files[@]}"
test "${#files[@]}" -eq 10

for f in "${files[@]}"; do
    test -s "$f"
    echo "PASS INPUT: $(basename "$f")"
done

fastqc \
  --threads 8 \
  --outdir "$OUTPUT" \
  "${files[@]}"

html_count=$(
    find "$OUTPUT" -maxdepth 1 -type f \
      -name '*_fastqc.html' | wc -l
)

zip_count=$(
    find "$OUTPUT" -maxdepth 1 -type f \
      -name '*_fastqc.zip' | wc -l
)

echo "FASTQC HTML REPORTS: $html_count"
echo "FASTQC ZIP REPORTS:  $zip_count"

test "$html_count" -eq 10
test "$zip_count" -eq 10

echo "UVC WGS POST-TRIM FASTQC: COMPLETE"
echo "FINISH: $(date)"
