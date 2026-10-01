#!/usr/bin/env bash
#SBATCH --job-name=GT_WGS_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --array=1-5%2
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/UVC/primary23/genotype_P23_%A_%a.out
#SBATCH --error=logs/UVC/primary23/genotype_P23_%A_%a.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"

CHR="${SLURM_ARRAY_TASK_ID}"
WORKSPACE="${ROOT}/results/wgs_variants/genomicsdb_primary23/chr${CHR}"
OUTDIR="${ROOT}/results/wgs_variants/joint_primary23"

OUT="${OUTDIR}/PRIMARY23_chr${CHR}.vcf.gz"
TBI="${OUT}.tbi"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load gatk/4.6.1.0

export TMPDIR="${SLURM_TMPDIR}"

echo "========================================"
echo "PRIMARY23 JOINT GENOTYPING"
echo "CHROMOSOME: ${CHR}"
echo "WORKSPACE:  ${WORKSPACE}"
echo "OUTPUT:     ${OUT}"
echo "HOST:       $(hostname)"
echo "START:      $(date)"
echo "TMPDIR:     ${TMPDIR}"
echo "========================================"

echo "===== INPUT VALIDATION ====="

for f in \
    "$REF" \
    "${REF}.fai" \
    "${ROOT}/reference/TAIR10/TAIR10_chr_all.dict" \
    "${WORKSPACE}/callset.json" \
    "${WORKSPACE}/vidmap.json" \
    "${WORKSPACE}/vcfheader.vcf"
do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty input: $f"
        exit 1
    fi

    printf "PASS\t%s\t%s bytes\n" \
        "$f" "$(stat -Lc '%s' "$f")"
done

if [[ -e "$OUT" || -e "$TBI" ]]; then
    echo "ERROR: output already exists for chromosome ${CHR}"
    ls -lh "$OUT" "$TBI" 2>/dev/null || true
    exit 1
fi

callset_count=$(
    grep -o '"sample_name"' "${WORKSPACE}/callset.json" |
    wc -l
)

echo "GENOMICSDB CALLSET COUNT: ${callset_count}"

if [[ "$callset_count" -ne 23 ]]; then
    echo "ERROR: expected 23 GenomicsDB callsets"
    exit 1
fi

echo "===== GENOTYPEGVCFS ====="

gatk \
    --java-options "-Xmx18g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://${WORKSPACE}" \
    -O "$OUT" \
    -L "$CHR" \
    --include-non-variant-sites false \
    --tmp-dir "$SLURM_TMPDIR"

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

sample_count=$(
    zgrep -m1 '^#CHROM' "$OUT" |
    awk -F '\t' '{print NF-9}'
)

variant_count=$(
    zgrep -vc '^#' "$OUT"
)

echo "VCF SAMPLE COUNT:  ${sample_count}"
echo "VCF VARIANT COUNT: ${variant_count}"

if [[ "$sample_count" -ne 23 ]]; then
    echo "ERROR: expected 23 samples in output VCF"
    exit 1
fi

if [[ "$variant_count" -eq 0 ]]; then
    echo "ERROR: output VCF contains no variant records"
    exit 1
fi

echo "===== VARIANT VALIDATION ====="

gatk \
    --java-options "-Xmx4g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    ValidateVariants \
    -R "$REF" \
    -V "$OUT" \
    -L "$CHR"

echo "VARIANT VALIDATION: PASS"
echo "PRIMARY23 JOINT GENOTYPING COMPLETE: CHR${CHR}"
echo "GENOMICSDB AND GVCF INPUTS RETAINED"
echo "FINISH: $(date)"
