#!/usr/bin/env bash
#SBATCH --job-name=QC_WGS_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/primary23/QC_raw_P23_%j.out
#SBATCH --error=logs/UVC/primary23/QC_raw_P23_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
VCF="${ROOT}/results/wgs_variants/joint_primary23/PRIMARY23_chr1-5.vcf.gz"
META="${ROOT}/metadata/wgs_primary23_groups.tsv"
OUTDIR="${ROOT}/results/wgs_variants/PRIMARY23_QC_raw"
PREFIX="${OUTDIR}/PRIMARY23_raw"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load bcftools/1.22
module load vcftools/0.1.16

echo "========================================"
echo "PRIMARY23 RAW JOINT-VCF QC"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

echo "===== SOFTWARE ====="
bcftools --version | head -1
vcftools --version

echo "===== INPUT VALIDATION ====="

test -s "$VCF"
test -s "${VCF}.tbi"
test -s "$META"
gzip -t "$VCF"

vcf_samples=$(bcftools query -l "$VCF" | wc -l)
metadata_samples=$(awk 'NR>1 {n++} END {print n+0}' "$META")
variants=$(bcftools view -H "$VCF" | wc -l)

echo "VCF SAMPLES:      $vcf_samples"
echo "METADATA SAMPLES: $metadata_samples"
echo "RAW VARIANTS:     $variants"

if [[ "$vcf_samples" -ne 23 ]]; then
    echo "ERROR: VCF does not contain 23 samples"
    exit 1
fi

if [[ "$metadata_samples" -ne 23 ]]; then
    echo "ERROR: metadata does not contain 23 samples"
    exit 1
fi

if [[ "$variants" -ne 51173 ]]; then
    echo "ERROR: expected 51173 variants"
    exit 1
fi

comm -3 \
  <(bcftools query -l "$VCF" | sort) \
  <(awk -F '\t' 'NR>1 {print $1}' "$META" | sort) \
  > "${PREFIX}.sample_disagreement.txt"

if [[ -s "${PREFIX}.sample_disagreement.txt" ]]; then
    echo "ERROR: VCF and metadata sample sets differ"
    cat "${PREFIX}.sample_disagreement.txt"
    exit 1
fi

echo "INPUT VALIDATION: PASS"

echo "===== HEADER INVENTORY ====="

bcftools view -h "$VCF" \
  > "${PREFIX}.header.vcf.txt"

bcftools query -l "$VCF" \
  > "${PREFIX}.samples.txt"

bcftools query -f '%CHROM\n' "$VCF" |
sort -V |
uniq -c \
  > "${PREFIX}.variants_by_chromosome.txt"

echo "===== BCFTOOLS GLOBAL AND SAMPLE STATISTICS ====="

bcftools stats \
  --threads 2 \
  --samples - \
  "$VCF" \
  > "${PREFIX}.bcftools_stats.txt"

grep '^SN' "${PREFIX}.bcftools_stats.txt" \
  > "${PREFIX}.summary_numbers.tsv" || true

grep '^PSC' "${PREFIX}.bcftools_stats.txt" \
  > "${PREFIX}.per_sample_counts.tsv" || true

grep '^TSTV' "${PREFIX}.bcftools_stats.txt" \
  > "${PREFIX}.transition_transversion.tsv" || true

grep '^IDD' "${PREFIX}.bcftools_stats.txt" \
  > "${PREFIX}.indel_distribution.tsv" || true

echo "===== VARIANT-TYPE COUNTS ====="

{
    printf "category\tcount\n"

    printf "all\t"
    bcftools view -H "$VCF" | wc -l

    printf "snps\t"
    bcftools view -v snps -H "$VCF" | wc -l

    printf "indels\t"
    bcftools view -v indels -H "$VCF" | wc -l

    printf "mnps\t"
    bcftools view -v mnps -H "$VCF" | wc -l

    printf "other\t"
    bcftools view -v other -H "$VCF" | wc -l

    printf "biallelic_snps\t"
    bcftools view -m2 -M2 -v snps -H "$VCF" | wc -l

    printf "multiallelic\t"
    bcftools view -m3 -H "$VCF" | wc -l
} > "${PREFIX}.variant_type_counts.tsv"

echo "===== VCFTOOLS MISSINGNESS ====="

vcftools \
  --gzvcf "$VCF" \
  --missing-indv \
  --out "${PREFIX}.missingness"

vcftools \
  --gzvcf "$VCF" \
  --missing-site \
  --out "${PREFIX}.site_missingness"

echo "===== VCFTOOLS DEPTH ====="

vcftools \
  --gzvcf "$VCF" \
  --depth \
  --out "${PREFIX}.depth"

vcftools \
  --gzvcf "$VCF" \
  --site-mean-depth \
  --out "${PREFIX}.site_depth"

echo "===== VCFTOOLS HETEROZYGOSITY ====="

vcftools \
  --gzvcf "$VCF" \
  --het \
  --out "${PREFIX}.heterozygosity"

echo "===== FORMAT-FIELD INVENTORY ====="

bcftools view -h "$VCF" |
grep '^##FORMAT=' \
  > "${PREFIX}.FORMAT_fields.txt"

echo "===== OUTPUT VALIDATION ====="

required=(
  "${PREFIX}.header.vcf.txt"
  "${PREFIX}.samples.txt"
  "${PREFIX}.variants_by_chromosome.txt"
  "${PREFIX}.bcftools_stats.txt"
  "${PREFIX}.summary_numbers.tsv"
  "${PREFIX}.per_sample_counts.tsv"
  "${PREFIX}.variant_type_counts.tsv"
  "${PREFIX}.missingness.imiss"
  "${PREFIX}.site_missingness.lmiss"
  "${PREFIX}.depth.idepth"
  "${PREFIX}.site_depth.ldepth.mean"
  "${PREFIX}.heterozygosity.het"
  "${PREFIX}.FORMAT_fields.txt"
)

for file in "${required[@]}"; do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: missing or empty QC file: $file"
        exit 1
    fi
done

echo "QC OUTPUT FILES: ${#required[@]}"
echo "OUTPUT VALIDATION: PASS"

echo
echo "===== VARIANT TYPES ====="
column -t "${PREFIX}.variant_type_counts.tsv"

echo
echo "===== PER-SAMPLE MISSINGNESS ====="
column -t "${PREFIX}.missingness.imiss"

echo
echo "===== PER-SAMPLE MEAN DEPTH ====="
column -t "${PREFIX}.depth.idepth"

echo
echo "===== PER-SAMPLE HETEROZYGOSITY ====="
column -t "${PREFIX}.heterozygosity.het"

echo
echo "PRIMARY23 RAW QC: COMPLETE"
echo "FINISH: $(date)"
