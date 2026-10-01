#!/usr/bin/env bash
#SBATCH --job-name=UVC_WGBS_trim
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-5%3
#SBATCH --output=logs/UVC/trimmed/UVC_trim_%A_%a.out
#SBATCH --error=logs/UVC/trimmed/UVC_trim_%A_%a.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
TRIM="${ROOT}/tools/trim_galore-2.3.0/trim_galore"

samples=(
    F25UV_7_WGBS
    F25UV_8_WGBS
    F25UV_10_WGBS
    F25UV_11_WGBS
    F25UV_13_WGBS
)

sample="${samples[$((SLURM_ARRAY_TASK_ID - 1))]}"

r1="${ROOT}/raw/UVC_normalized/${sample}_1.fastq.gz"
r2="${ROOT}/raw/UVC_normalized/${sample}_2.fastq.gz"

outdir="${ROOT}/trimmed/wgbs_all"
out1="${outdir}/${sample}_1_val_1.fq.gz"
out2="${outdir}/${sample}_2_val_2.fq.gz"
trimlog="${ROOT}/logs/UVC/trimmed/${sample}.trim_galore.log"

echo "========================================"
echo "SAMPLE: ${sample}"
echo "START:  $(date)"
echo "HOST:   $(hostname)"
echo "========================================"

[[ -x "$TRIM" ]] || {
    echo "ERROR: Trim Galore executable missing: $TRIM"
    exit 1
}

for f in "$r1" "$r2"; do
    [[ -s "$f" ]] || {
        echo "ERROR: input missing or empty: $f"
        exit 1
    }
done

mkdir -p "$outdir"

# Safely accept an already completed sample.
if [[ -s "$out1" && -s "$out2" ]] &&
   gzip -t "$out1" &&
   gzip -t "$out2"; then
    echo "ALREADY COMPLETE AND VALID: ${sample}"
    exit 0
fi

# Do not overwrite ambiguous partial outputs automatically.
if [[ -e "$out1" || -e "$out2" ]]; then
    echo "ERROR: partial or invalid output already exists for ${sample}"
    echo "Inspect before rerunning:"
    echo "$out1"
    echo "$out2"
    exit 1
fi

"$TRIM" \
    --paired \
    --quality 20 \
    --length 20 \
    --cores 2 \
    --output_dir "$outdir" \
    "$r1" "$r2" \
    > "$trimlog" 2>&1

for f in "$out1" "$out2"; do
    [[ -s "$f" ]] || {
        echo "ERROR: expected trimmed output missing: $f"
        exit 1
    }

    gzip -t "$f" || {
        echo "ERROR: gzip integrity failure: $f"
        exit 1
    }
done

echo "TRIMMED OUTPUTS VALIDATED"
echo "R1: $(stat -c '%s bytes' "$out1")"
echo "R2: $(stat -c '%s bytes' "$out2")"

echo "FINISH: $(date)"
echo "COMPLETE: ${sample}"
