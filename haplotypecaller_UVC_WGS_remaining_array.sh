#!/usr/bin/env bash
#SBATCH --job-name=UVC_HC_prod
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --array=1-4%2
#SBATCH --output=logs/UVC/haplotypecaller/UVC_HC_prod_%A_%a.out
#SBATCH --error=logs/UVC/haplotypecaller/UVC_HC_prod_%A_%a.err

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
BAM="${ROOT}/results/wgs_alignment/UVC/${SAMPLE}.markdup.bam"
BAI="${BAM}.bai"

OUTDIR="${ROOT}/results/wgs_variants/UVC_gvcf"
OUT="${OUTDIR}/${SAMPLE}.g.vcf.gz"
TBI="${OUT}.tbi"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load gatk/4.6.1.0
module load samtools/1.22.1

export TMPDIR="${SLURM_TMPDIR}"

echo "========================================"
echo "UVC WGS HAPLOTYPECALLER PRODUCTION"
echo "SAMPLE: ${SAMPLE}"
echo "TASK:   ${SLURM_ARRAY_TASK_ID}"
echo "START:  $(date)"
echo "HOST:   $(hostname)"
echo "TMPDIR: ${TMPDIR}"
echo "========================================"

echo "===== SOFTWARE ====="
gatk --version
samtools --version | head -1

echo "===== INPUT VALIDATION ====="

for f in \
    "$REF" \
    "${REF}.fai" \
    "${ROOT}/reference/TAIR10/TAIR10_chr_all.dict" \
    "$BAM" \
    "$BAI"
do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty input: $f"
        exit 1
    fi

    printf "PASS\t%s\t%s bytes\n" \
        "$f" "$(stat -Lc '%s' "$f")"
done

if [[ -e "$OUT" || -e "$TBI" ]]; then
    echo "ERROR: output already exists for ${SAMPLE}"
    ls -lh "$OUT" "$TBI" 2>/dev/null || true
    exit 1
fi

samtools quickcheck -v "$BAM"

observed_sample=$(
    samtools view -H "$BAM" |
    awk -F '\t' '
        $1=="@RG" {
            for (i=2; i<=NF; i++) {
                if ($i ~ /^SM:/) {
                    sub(/^SM:/, "", $i)
                    print $i
                }
            }
        }
    ' |
    sort -u
)

if [[ "$observed_sample" != "$SAMPLE" ]]; then
    echo "ERROR: BAM sample is '${observed_sample}', expected '${SAMPLE}'"
    exit 1
fi

echo "BAM SAMPLE NAME: ${observed_sample}"
echo "INPUT VALIDATION: PASS"

echo "===== HAPLOTYPECALLER ====="

gatk \
    --java-options "-Xmx18g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    HaplotypeCaller \
    -R "$REF" \
    -I "$BAM" \
    -O "$OUT" \
    -ERC GVCF \
    --native-pair-hmm-threads 8 \
    --tmp-dir "$SLURM_TMPDIR" \
    -L 1 \
    -L 2 \
    -L 3 \
    -L 4 \
    -L 5

echo "===== OUTPUT VALIDATION ====="

for f in "$OUT" "$TBI"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty output: $f"
        exit 1
    fi

    printf "PASS\t%s\t%s bytes\n" \
        "$f" "$(stat -Lc '%s' "$f")"
done

gzip -t "$OUT"
echo "GZIP VALIDATION: PASS"

header_sample=$(
    zgrep -m1 '^#CHROM' "$OUT" |
    awk -F '\t' '{print $10}'
)

if [[ "$header_sample" != "$SAMPLE" ]]; then
    echo "ERROR: gVCF sample is '${header_sample}', expected '${SAMPLE}'"
    exit 1
fi

echo "GVCF SAMPLE NAME: ${header_sample}"

echo "===== NUCLEAR GVCF VALIDATION ====="

gatk \
    --java-options "-Xmx4g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    ValidateVariants \
    -R "$REF" \
    -V "$OUT" \
    -gvcf \
    -L 1 \
    -L 2 \
    -L 3 \
    -L 4 \
    -L 5

echo "NUCLEAR GVCF VALIDATION: PASS"

echo "===== OUTPUT SUMMARY ====="
printf "GVCF:\t%s bytes\n" "$(stat -Lc '%s' "$OUT")"
printf "INDEX:\t%s bytes\n" "$(stat -Lc '%s' "$TBI")"

echo "BAM AND FASTQ INPUTS RETAINED"
echo "UVC HAPLOTYPECALLER COMPLETE: ${SAMPLE}"
echo "FINISH: $(date)"
