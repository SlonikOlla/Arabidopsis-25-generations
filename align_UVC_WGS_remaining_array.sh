#!/usr/bin/env bash
#SBATCH --job-name=UVC_WGS_prod
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-4%2
#SBATCH --output=logs/UVC/wgs_alignment/UVC_WGS_prod_%A_%a.out
#SBATCH --error=logs/UVC/wgs_alignment/UVC_WGS_prod_%A_%a.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"

samples=(
    F25UV_8_WGS
    F25UV_10_WGS
    F25UV_11_WGS
    F25UV_13_WGS
)

SAMPLE="${samples[$((SLURM_ARRAY_TASK_ID - 1))]}"

REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"

R1="${ROOT}/trimmed/UVC_WGS/${SAMPLE}_1_val_1.fq.gz"
R2="${ROOT}/trimmed/UVC_WGS/${SAMPLE}_2_val_2.fq.gz"

OUTDIR="${ROOT}/results/wgs_alignment/UVC"
QCDIR="${OUTDIR}/qc"
METRICDIR="${OUTDIR}/metrics"

OUT="${OUTDIR}/${SAMPLE}.markdup.bam"
BAI="${OUT}.bai"
FLAGSTAT="${QCDIR}/${SAMPLE}.flagstat.txt"
IDXSTATS="${QCDIR}/${SAMPLE}.idxstats.txt"
METRICS="${METRICDIR}/${SAMPLE}.duplicate_metrics.txt"

SORTED_TMP="${SLURM_TMPDIR}/${SAMPLE}.sorted.bam"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load bwa-mem2/2.2.1
module load samtools/1.22.1
module load picard/3.1.0

mkdir -p "$OUTDIR" "$QCDIR" "$METRICDIR"

echo "========================================"
echo "UVC WGS PRODUCTION"
echo "SAMPLE: $SAMPLE"
echo "TASK:   $SLURM_ARRAY_TASK_ID"
echo "START:  $(date)"
echo "HOST:   $(hostname)"
echo "========================================"

echo "===== SOFTWARE ====="
bwa-mem2 version 2>&1 | head -3 || true
samtools --version | head -1

echo "===== INPUT VALIDATION ====="

for f in "$REF" "$R1" "$R2"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty input: $f"
        exit 1
    fi

    printf "PASS\t%s\t%s bytes\n" \
        "$f" "$(stat -Lc '%s' "$f")"
done

if [[ -e "$OUT" || -e "$BAI" ]]; then
    echo "ERROR: output already exists for $SAMPLE"
    echo "$OUT"
    echo "$BAI"
    exit 1
fi

echo "===== ALIGNMENT AND SORTING ====="

bwa-mem2 mem \
    -t 8 \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:UVC_WGS\tPU:${SAMPLE}" \
    "$REF" "$R1" "$R2" |
samtools sort \
    -@ 4 \
    -m 3G \
    -T "${SLURM_TMPDIR}/${SAMPLE}.sort" \
    -o "$SORTED_TMP" -

samtools quickcheck -v "$SORTED_TMP"

echo "SORTED BAM VALIDATED"

echo "===== DUPLICATE MARKING ====="

java -Xmx24g \
    -jar "$EBROOTPICARD/picard.jar" \
    MarkDuplicates \
    INPUT="$SORTED_TMP" \
    OUTPUT="$OUT" \
    METRICS_FILE="$METRICS" \
    TMP_DIR="$SLURM_TMPDIR" \
    CREATE_INDEX=false \
    VALIDATION_STRINGENCY=SILENT \
    ASSUME_SORT_ORDER=coordinate

test -s "$OUT"
test -s "$METRICS"

echo "===== INDEX AND QC ====="

samtools index -@ 4 "$OUT"
samtools quickcheck -v "$OUT"
samtools flagstat -@ 4 "$OUT" > "$FLAGSTAT"
samtools idxstats "$OUT" > "$IDXSTATS"

test -s "$BAI"
test -s "$FLAGSTAT"
test -s "$IDXSTATS"

echo "===== OUTPUT VALIDATION ====="
printf "BAM:\t%s bytes\n" "$(stat -c '%s' "$OUT")"
printf "BAI:\t%s bytes\n" "$(stat -c '%s' "$BAI")"
printf "METRICS:\t%s bytes\n" "$(stat -c '%s' "$METRICS")"

echo "===== ALIGNMENT SUMMARY ====="
cat "$FLAGSTAT"

echo "===== DUPLICATION SUMMARY ====="
awk -F '\t' '
    /^LIBRARY/ {
        getline
        print "LIBRARY\tREAD_PAIRS_EXAMINED\tREAD_PAIR_DUPLICATES\tPERCENT_DUPLICATION"
        print $1 "\t" $3 "\t" $7 "\t" $9
        exit
    }
' "$METRICS"

echo
echo "RAW AND TRIMMED FASTQS RETAINED"
echo "COMPLETE: $SAMPLE"
echo "FINISH: $(date)"
