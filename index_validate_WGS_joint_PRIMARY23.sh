#!/usr/bin/env bash
#SBATCH --job-name=Index_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/UVC/primary23/index_P23_%j.out
#SBATCH --error=logs/UVC/primary23/index_P23_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"
VCF="${ROOT}/results/wgs_variants/joint_primary23/PRIMARY23_chr1-5.vcf.gz"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load gatk/4.6.1.0

export JAVA_TOOL_OPTIONS="-Xmx3g"

echo "========================================"
echo "PRIMARY23 INDEX AND VALIDATION"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

echo "===== VCF PREFLIGHT ====="

test -s "$VCF"
gzip -t "$VCF"

sample_count=$(
    zgrep -m1 '^#CHROM' "$VCF" |
    awk -F '\t' '{print NF-9}'
)

variant_count=$(zgrep -vc '^#' "$VCF")

echo "VCF SIZE:          $(stat -Lc '%s' "$VCF") bytes"
echo "VCF SAMPLE COUNT:  $sample_count"
echo "VCF VARIANT COUNT: $variant_count"

if [[ "$sample_count" -ne 23 ]]; then
    echo "ERROR: expected 23 samples"
    exit 1
fi

if [[ "$variant_count" -ne 51173 ]]; then
    echo "ERROR: expected 51173 variants"
    exit 1
fi

if [[ -e "${VCF}.tbi" ]]; then
    echo "ERROR: index already exists: ${VCF}.tbi"
    exit 1
fi

echo "VCF PREFLIGHT: PASS"

echo "===== CREATING INDEX ====="

gatk --java-options \
    "-Xmx3g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    IndexFeatureFile \
    -I "$VCF"

test -s "${VCF}.tbi"

echo "INDEX CREATION: PASS"
echo "INDEX SIZE: $(stat -Lc '%s' "${VCF}.tbi") bytes"

echo "===== GATK VALIDATION ====="

gatk --java-options \
    "-Xmx3g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    ValidateVariants \
    -R "$REF" \
    -V "$VCF" \
    -L 1 \
    -L 2 \
    -L 3 \
    -L 4 \
    -L 5

echo "GZIP VALIDATION: PASS"
echo "INDEX VALIDATION: PASS"
echo "VARIANT VALIDATION: PASS"
echo "PRIMARY23 COMBINED VCF: COMPLETE"
echo "OUTPUT: $VCF"
echo "FINISH: $(date)"
