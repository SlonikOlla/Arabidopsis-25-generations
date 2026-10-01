#!/usr/bin/env bash
#SBATCH --job-name=Gather_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/primary23/gather_P23_%j.out
#SBATCH --error=logs/UVC/primary23/gather_P23_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"
INDIR="${ROOT}/results/wgs_variants/joint_primary23"
OUT="${INDIR}/PRIMARY23_chr1-5.vcf.gz"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load gatk/4.6.1.0

export JAVA_TOOL_OPTIONS="-Xmx6g"

echo "========================================"
echo "PRIMARY23 VCF GATHER"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

echo "===== INPUT VALIDATION ====="

for chr in 1 2 3 4 5; do
    vcf="${INDIR}/PRIMARY23_chr${chr}.vcf.gz"

    test -s "$vcf"
    test -s "${vcf}.tbi"
    gzip -t "$vcf"

    samples=$(
        zgrep -m1 '^#CHROM' "$vcf" |
        awk -F '\t' '{print NF-9}'
    )

    if [[ "$samples" -ne 23 ]]; then
        echo "ERROR: chromosome ${chr} contains ${samples} samples"
        exit 1
    fi

    echo "PASS: CHR${chr}, ${samples} samples"
done

echo "===== SAMPLE-ORDER VALIDATION ====="

reference_header=$(
    zgrep -m1 '^#CHROM' "${INDIR}/PRIMARY23_chr1.vcf.gz"
)

for chr in 2 3 4 5; do
    current_header=$(
        zgrep -m1 '^#CHROM' "${INDIR}/PRIMARY23_chr${chr}.vcf.gz"
    )

    if [[ "$current_header" != "$reference_header" ]]; then
        echo "ERROR: sample order differs on chromosome ${chr}"
        exit 1
    fi
done

echo "SAMPLE ORDER: PASS"

if [[ -e "$OUT" || -e "${OUT}.tbi" ]]; then
    echo "ERROR: output already exists: $OUT"
    exit 1
fi

echo "===== GATHERING CHROMOSOMES ====="

gatk --java-options "-Xmx6g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    GatherVcfs \
    -I "${INDIR}/PRIMARY23_chr1.vcf.gz" \
    -I "${INDIR}/PRIMARY23_chr2.vcf.gz" \
    -I "${INDIR}/PRIMARY23_chr3.vcf.gz" \
    -I "${INDIR}/PRIMARY23_chr4.vcf.gz" \
    -I "${INDIR}/PRIMARY23_chr5.vcf.gz" \
    -O "$OUT" \
    --CREATE_INDEX true

echo "===== OUTPUT VALIDATION ====="

test -s "$OUT"
test -s "${OUT}.tbi"
gzip -t "$OUT"

sample_count=$(
    zgrep -m1 '^#CHROM' "$OUT" |
    awk -F '\t' '{print NF-9}'
)

variant_count=$(zgrep -vc '^#' "$OUT")

if [[ "$sample_count" -ne 23 ]]; then
    echo "ERROR: gathered VCF contains ${sample_count} samples"
    exit 1
fi

if [[ "$variant_count" -ne 51173 ]]; then
    echo "ERROR: expected 51173 variants but found ${variant_count}"
    exit 1
fi

gatk --java-options "-Xmx6g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    ValidateVariants \
    -R "$REF" \
    -V "$OUT" \
    -L 1 \
    -L 2 \
    -L 3 \
    -L 4 \
    -L 5

echo "GZIP VALIDATION: PASS"
echo "VCF SAMPLE COUNT:  $sample_count"
echo "VCF VARIANT COUNT: $variant_count"
echo "VARIANT VALIDATION: PASS"
echo "PRIMARY23 GATHER COMPLETE"
echo "OUTPUT: $OUT"
echo "FINISH: $(date)"
