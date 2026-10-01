#!/usr/bin/env bash
#SBATCH --job-name=PCA5kb_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/primary23/PCA_5kb_P23_%j.out
#SBATCH --error=logs/UVC/primary23/PCA_5kb_P23_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
PCADIR="${ROOT}/results/wgs_variants/PRIMARY23_PCA"
FILTERDIR="${ROOT}/results/wgs_variants/PRIMARY23_filtered"

PCA_VCF="${FILTERDIR}/PRIMARY23_PCA_SNPs.vcf.gz"
PLINK_PREFIX="${PCADIR}/PRIMARY23_PCA_input"

THIN_LIST="${PCADIR}/PRIMARY23_PCA_5kb_SNPs.txt"
PCA_ALL="${PCADIR}/PRIMARY23_all23_5kb"
PCA_NO65="${PCADIR}/PRIMARY23_without_SRR17867665_5kb"
REMOVE65="${PCADIR}/remove_SRR17867665_5kb.txt"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load bcftools/1.22
module load plink/2.0.0-a.6.32

echo "========================================"
echo "PRIMARY23 PCA RECOVERY: 5-KB THINNING"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

echo "===== INPUT VALIDATION ====="

test -s "$PCA_VCF"
test -s "${PCA_VCF}.tbi"
test -s "${PLINK_PREFIX}.pgen"
test -s "${PLINK_PREFIX}.pvar"
test -s "${PLINK_PREFIX}.psam"

gzip -t "$PCA_VCF"

vcf_samples=$(bcftools query -l "$PCA_VCF" | wc -l)
vcf_variants=$(bcftools view -H "$PCA_VCF" | wc -l)
plink_samples=$(awk 'NR>1 {n++} END {print n+0}' "${PLINK_PREFIX}.psam")
plink_variants=$(grep -vc '^#' "${PLINK_PREFIX}.pvar")

echo "PCA VCF SAMPLES:    $vcf_samples"
echo "PCA VCF VARIANTS:   $vcf_variants"
echo "PLINK SAMPLES:      $plink_samples"
echo "PLINK VARIANTS:     $plink_variants"

if [[ "$vcf_samples" -ne 23 ]]; then
    echo "ERROR: expected 23 VCF samples"
    exit 1
fi

if [[ "$plink_samples" -ne 23 ]]; then
    echo "ERROR: expected 23 PLINK samples"
    exit 1
fi

if [[ "$vcf_variants" -ne "$plink_variants" ]]; then
    echo "ERROR: VCF and PLINK variant counts differ"
    exit 1
fi

echo "INPUT VALIDATION: PASS"

echo "===== DETERMINISTIC 5-KB PHYSICAL THINNING ====="

awk '
BEGIN {
    chromosome = ""
    last_position = -1000000000
}

/^#/ {
    next
}

{
    current_chromosome = $1
    current_position = $2
    variant_id = $3

    if (
        current_chromosome != chromosome ||
        current_position - last_position >= 5000
    ) {
        print variant_id
        chromosome = current_chromosome
        last_position = current_position
    }
}
' "${PLINK_PREFIX}.pvar" > "$THIN_LIST"

test -s "$THIN_LIST"

thin_count=$(wc -l < "$THIN_LIST")

echo "INPUT PCA SNPS:     $plink_variants"
echo "5-KB THINNED SNPS:  $thin_count"

if [[ "$thin_count" -lt 100 ]]; then
    echo "ERROR: fewer than 100 physically thinned SNPs"
    exit 1
fi

echo "===== PCA: ALL 23 SAMPLES ====="

plink2 \
  --pfile "$PLINK_PREFIX" \
  --allow-extra-chr \
  --extract "$THIN_LIST" \
  --pca 10 allele-wts \
  --out "$PCA_ALL"

test -s "${PCA_ALL}.eigenvec"
test -s "${PCA_ALL}.eigenval"
test -s "${PCA_ALL}.eigenvec.allele"

all_samples=$(
    awk 'NR>1 {n++} END {print n+0}' \
      "${PCA_ALL}.eigenvec"
)

echo "ALL-SAMPLE PCA RECORDS: $all_samples"

if [[ "$all_samples" -ne 23 ]]; then
    echo "ERROR: expected 23 samples in primary PCA"
    exit 1
fi

echo "===== PCA WITHOUT SRR17867665 ====="

printf "SRR17867665\tSRR17867665\n" > "$REMOVE65"

plink2 \
  --pfile "$PLINK_PREFIX" \
  --allow-extra-chr \
  --extract "$THIN_LIST" \
  --remove "$REMOVE65" \
  --pca 10 allele-wts \
  --out "$PCA_NO65"

test -s "${PCA_NO65}.eigenvec"
test -s "${PCA_NO65}.eigenval"
test -s "${PCA_NO65}.eigenvec.allele"

no65_samples=$(
    awk 'NR>1 {n++} END {print n+0}' \
      "${PCA_NO65}.eigenvec"
)

echo "NO-SRR17867665 PCA RECORDS: $no65_samples"

if [[ "$no65_samples" -ne 22 ]]; then
    echo "ERROR: expected 22 samples in sensitivity PCA"
    exit 1
fi

echo "===== EIGENVALUES: ALL 23 ====="
nl -ba "${PCA_ALL}.eigenval"

echo "===== EIGENVALUES: WITHOUT SRR17867665 ====="
nl -ba "${PCA_NO65}.eigenval"

echo "===== FINAL VALIDATION ====="
echo "FILTERED PCA-ELIGIBLE SNPS: $plink_variants"
echo "5-KB THINNED PCA SNPS:      $thin_count"
echo "PRIMARY PCA SAMPLES:         $all_samples"
echo "SENSITIVITY PCA SAMPLES:     $no65_samples"
echo "SRR17867665 RETAINED IN PRIMARY PCA"
echo "NO PLINK --bad-ld OVERRIDE USED"
echo "PRIMARY23 5-KB PCA: COMPLETE"
echo "FINISH: $(date)"
